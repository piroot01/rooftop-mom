classdef ComputeMagPartChunked < mom.ImpedanceMatrix.ComputeStrategy
    properties (Access = private)
        computeContext_ mom.ImpedanceMatrix.ComputeChunkedContext;
    end

    methods
        function obj = ComputeMagPartChunked(computeContext)
            obj.computeContext_ = computeContext;
        end
    end

    methods
        function Z = compute(obj, strategyContext, varargin)
            k = strategyContext.constants.k;
            Z0 = strategyContext.constants.Z0;
            polaritiesCount = strategyContext.constants.polaritiesCount;

            order = strategyContext.quadrature.getOrder();
            basisFunctionCount = strategyContext.basisFunctions.getCount();

            if nargin == 3
                chunks = varargin{1};
            else
                chunks = obj.computeContext_.generateChunksAllInteractions(order, basisFunctionCount);
            end

            chunkCount = numel(chunks);
            blockSize = order^2;

            quadratureWeights = strategyContext.quadrature.getWeights();
            quadraturePoints = strategyContext.rooftopBasisFunctions.getQuadraturePoints();
            rooftopPoints = strategyContext.rooftopBasisFunctions.getRooftopPoints();

            slices = obj.computeContext_.generateSlicesFromChunks(chunks, blockSize);

            singularityMask = obj.getSingularityMask(strategyContext);
            singularityMaskNonZeroDot = false(size(singularityMask));
            IDMask = obj.getIDMask(strategyContext);
            jacobians = obj.getJacobians(strategyContext);

            Z = zeros(basisFunctionCount, basisFunctionCount, polaritiesCount);
            for n = 1:chunkCount
                % for debugging
                fprintf("Chunk ID=%d\n", n);

                chunk = chunks{n};
                slice = slices{n};

                rowSlice = slice.rowIndicesStart:slice.rowIndicesEnd;
                colSlice = slice.columnIndicesStart:slice.columnIndicesEnd;

                chunkContext.chunkSize = [numel(chunk.rowIndices), numel(chunk.columnIndices)];
                chunkContext.rowQuadraturePoints = quadraturePoints(:, rowSlice, :);
                chunkContext.columnQuadraturePoints = quadraturePoints(:, colSlice, :);
                chunkContext.rowRooftopPoints = rooftopPoints(:, rowSlice, :);
                chunkContext.columnRooftopPoints = rooftopPoints(:, colSlice, :);

                chunkContext.singularityMaskChunk = singularityMask(chunk.rowIndices, chunk.columnIndices, :);

                [ZBlock, singularityMaskChunkNonZeroDot] = obj.computeZBlock(strategyContext, chunkContext, quadratureWeights, blockSize);
                Z(chunk.rowIndices, chunk.columnIndices, 1:polaritiesCount) = ZBlock;
                singularityMaskNonZeroDot(chunk.rowIndices, chunk.columnIndices, :) = singularityMaskChunkNonZeroDot;
            end

            Z = Z + obj.computeMagSingularities(strategyContext, singularityMaskNonZeroDot, IDMask);

            Z = 1i * k * Z0 / (4 * pi) .* Z .* jacobians;

            Z = sum(Z, 3);
        end
    end

    methods (Access = private)
        function [ZBlock, singularityMaskChunkNonZeroDot] = computeZBlock(~, strategyContext, chunkContext, quadratureWeights, blockSize)
            k = strategyContext.constants.k;
            weights = repmat(quadratureWeights, chunkContext.chunkSize);

            p1 = strategyContext.constants.polarities(:, 1);
            p2 = strategyContext.constants.polarities(:, 2);

            rowQuadraturePointsPolarity = chunkContext.rowQuadraturePoints(:, :, p1);
            columnQuadraturePointsPolarity = chunkContext.columnQuadraturePoints(:, :, p2);
            rowRooftopPointsPolarity = chunkContext.rowRooftopPoints(:, :, p1);
            columnRooftopPointsPolarity = chunkContext.columnRooftopPoints(:, :, p2);

            bigSingularityMask = repelem(chunkContext.singularityMaskChunk, blockSize, blockSize);

            distanceMatrix = mom.DistanceMatrix.compute(rowQuadraturePointsPolarity, columnQuadraturePointsPolarity);
            distanceMatrixZero = (distanceMatrix < eps);

            rooftopPointsDot = pagemtimes(permute(rowRooftopPointsPolarity, [2, 1, 3]), columnRooftopPointsPolarity);
            rooftopPointsDotZero = (rooftopPointsDot < eps);

            bigSingularityMaskNonZeroDot = ~rooftopPointsDotZero & bigSingularityMask;
            singularityMaskChunkNonZeroDot = bigSingularityMaskNonZeroDot(1:blockSize:end, 1:blockSize:end, :);

            greensFunction = exp(-1i * k * distanceMatrix) ./ distanceMatrix;
            greensFunction(bigSingularityMaskNonZeroDot) = greensFunction(bigSingularityMaskNonZeroDot) - 1 ./ distanceMatrix(bigSingularityMaskNonZeroDot);

            integrands = rooftopPointsDot .* greensFunction;
            integrands(distanceMatrixZero) = rooftopPointsDot(distanceMatrixZero) * (-1i) * k;

            integrals = integrands .* weights;
            integralsReshaped = reshape(integrals, blockSize, chunkContext.chunkSize(1), blockSize, chunkContext.chunkSize(2), strategyContext.constants.polaritiesCount);
            integralsSummed = sum(sum(integralsReshaped, 1), 3);
            ZBlock = reshape(integralsSummed, [chunkContext.chunkSize, strategyContext.constants.polaritiesCount]);
        end

        function singularities = computeMagSingularities(~, strategyContext, singularityMask, IDMask)
            [rooftopSizes, rooftopSizeIndices] = strategyContext.rooftopBasisFunctions.getUniqueRooftopSizes();
            [rowCount, colCount, ~] = size(singularityMask);
            rooftopSizeIndicesMatrix = repmat(rooftopSizeIndices, 1, colCount);

            singularityValuesType1 = mom.ImpedanceMatrix.singularity.type1(rooftopSizes(:,1), rooftopSizes(:,2));
            singularityValuesType2 = mom.ImpedanceMatrix.singularity.type2(rooftopSizes(:,1), rooftopSizes(:,2));

            singularityValuesType1Matrix = singularityValuesType1(rooftopSizeIndicesMatrix);
            singularityValuesType2Matrix = singularityValuesType2(rooftopSizeIndicesMatrix);

            validIndices = (IDMask > 0);
            combinedMask = bsxfun(@and, singularityMask, validIndices);
            singularities = zeros(rowCount, colCount, strategyContext.constants.polaritiesCount);

            samePolarities = [1, 4];
            differentPolarities = [2, 3];

            singularities(:, :, samePolarities) = combinedMask(:, :, samePolarities) .* repmat(singularityValuesType1Matrix, [1, 1, 2]);
            singularities(:, :, differentPolarities) = combinedMask(:, :, differentPolarities) .* repmat(singularityValuesType2Matrix, [1, 1, 2]);
        end
    end
end
