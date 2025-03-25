classdef ComputeElPartChunked < mom.ImpedanceMatrix.ComputeStrategy
    properties (Access = private)
        computeContext_ mom.ImpedanceMatrix.ComputeChunkedContext;
    end

    methods
        function obj = ComputeElPartChunked(computeContext)
            obj.computeContext_ = computeContext;
        end
    end

    methods
        function Z = compute(obj, strategyContext, varargin)
            k = strategyContext.constants.k;
            Z0 = strategyContext.constants.Z0;
            polaritiesCount = strategyContext.constants.polaritiesCount;
            p1 = strategyContext.constants.polarities(:, 1);
            p2 = strategyContext.constants.polarities(:, 2);

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
            rooftopDivPointsDot = mom.ImpedanceMatrix.ComputeElPartChunked.computeRooftopDivPointsDot(strategyContext);

            slices = obj.computeContext_.generateSlicesFromChunks(chunks, blockSize);

            singularityMask = obj.getSingularityMask(strategyContext);
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
                chunkContext.rooftopDivPointsDot = rooftopDivPointsDot;
                chunkContext.singularityMaskChunk = singularityMask(chunk.rowIndices, chunk.columnIndices, :);

                rowQuadraturePointsPolarity = quadraturePoints(:, rowSlice, p1);
                columnQuadraturePointsPolarity = quadraturePoints(:, colSlice, p2);

                chunkContext.integralKernel = mom.ImpedanceMatrix.IntegralKernel.compute(strategyContext, rowQuadraturePointsPolarity, columnQuadraturePointsPolarity);

                ZBlock = mom.ImpedanceMatrix.ComputeElPartChunked.computeZBlock(strategyContext, chunkContext, quadratureWeights, blockSize);
                Z(chunk.rowIndices, chunk.columnIndices, 1:polaritiesCount) = ZBlock;
            end

            Z = (Z + mom.ImpedanceMatrix.ComputeElPartChunked.computeElSingularities(strategyContext, singularityMask, IDMask)) .* rooftopDivPointsDot;

            Z = 1i * Z0 / (4 * pi * k) .* Z .* jacobians;

            Z = sum(Z, 3);
        end
    end

    methods (Access = public, Static)
        function ZBlock = computeZBlock(strategyContext, chunkContext, quadratureWeights, blockSize)
            k = strategyContext.constants.k;
            weights = repmat(quadratureWeights, chunkContext.chunkSize);

            bigSingularityMask = repelem(chunkContext.singularityMaskChunk, blockSize, blockSize);

            distanceMatrix = chunkContext.integralKernel.distanceMatrix;
            distanceMatrixZero = (distanceMatrix < eps);

            greensFunction = chunkContext.integralKernel.greensFunction;
            greensFunction(bigSingularityMask) = greensFunction(bigSingularityMask) - 1 ./ distanceMatrix(bigSingularityMask);

            integrands = greensFunction;
            integrands(distanceMatrixZero) = (-1i) * k;

            integrals = integrands .* weights;
            integralsReshaped = reshape(integrals, blockSize, chunkContext.chunkSize(1), blockSize, chunkContext.chunkSize(2), strategyContext.constants.polaritiesCount);
            integralsSummed = sum(sum(integralsReshaped, 1), 3);
            ZBlock = reshape(integralsSummed, [chunkContext.chunkSize, strategyContext.constants.polaritiesCount]);
        end

        function singularities = computeElSingularities(strategyContext, singularityMask, IDMask)
            [rooftopSizes, rooftopSizeIndices] = strategyContext.rooftopBasisFunctions.getUniqueRooftopSizes();
            [~, colCount, ~] = size(singularityMask);
            rooftopSizeIndicesMatrix = repmat(rooftopSizeIndices, 1, colCount);

            singularityValuesType3 = mom.ImpedanceMatrix.singularity.type3(rooftopSizes(:,1), rooftopSizes(:,2));

            singularityValuesType3Matrix = singularityValuesType3(rooftopSizeIndicesMatrix);

            validIndices = (IDMask > 0);
            combinedMask = bsxfun(@and, singularityMask, validIndices);

            singularities = combinedMask .* repmat(singularityValuesType3Matrix, [1, 1, strategyContext.constants.polaritiesCount]);
        end

        function result = computeRooftopDivPointsDot(strategyContext)
            rooftopDivPoints = strategyContext.rooftopBasisFunctions.getRooftopDivPoints();

            p1 = strategyContext.constants.polarities(:, 1);
            p2 = strategyContext.constants.polarities(:, 2);

            rowRooftopDivPointsPolarity = rooftopDivPoints(:, :, p1);
            columnRooftopDivPointsPolarity = rooftopDivPoints(:, :, p2);

            result = pagemtimes(permute(rowRooftopDivPointsPolarity, [2, 1, 3]), columnRooftopDivPointsPolarity);
        end
    end
end
