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

            Z = zeros(basisFunctionCount);
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
                Z(chunk.rowIndices, chunk.columnIndices) = ZBlock;
                singularityMaskNonZeroDot(chunk.rowIndices, chunk.columnIndices, :) = singularityMaskChunkNonZeroDot;
            end

            Z = Z + obj.computeMagSingularities(strategyContext, singularityMaskNonZeroDot, IDMask);

            Z = 1i * k * Z0 / (4 * pi) .* Z .* jacobians;
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
            integralsSummedReshaped = reshape(integralsSummed, [chunkContext.chunkSize, strategyContext.constants.polaritiesCount]);

            ZBlock = sum(integralsSummedReshaped, 3);
        end

        function singularities = computeMagSingularities(obj, strategyContext, singularityMask, IDMask)
            IDCount = strategyContext.basisFunctions.getIDCount();
            [rooftopSizes, rooftopSizeIndices] = strategyContext.rooftopBasisFunctions.getUniqueRooftopSizes();
            singularities = zeros(size(singularityMask));

            for n = 1:strategyContext.constants.polaritiesCount
                singularitiesPolarity = singularities(:, :, n);
                singularityValues = obj.evalMagSingularities(rooftopSizes, n);

                for id = 1:IDCount
                    matchIDMask = (id == IDMask);
                    matchIDMaskAndSingularityMask = matchIDMask & singularityMask(:, :, n);
                    [indices, ~] = find(matchIDMaskAndSingularityMask);
                    singularityValueMask = rooftopSizeIndices(indices);
                    singularitiesPolarity(matchIDMaskAndSingularityMask) = singularitiesPolarity(matchIDMaskAndSingularityMask) + singularityValues(singularityValueMask);
                end

                singularities(:, :, n) = singularitiesPolarity;
            end

            singularities = sum(singularities, 3);
        end

        function singularities = evalMagSingularities(~, xs, polarity)
            if (polarity == 1 || polarity == 4)
                singularities = mom.ImpedanceMatrix.singularity.type1(xs(:, 1), xs(:, 2));
            else
                singularities = mom.ImpedanceMatrix.singularity.type2(xs(:, 1), xs(:, 2));
            end
        end
    end
end
