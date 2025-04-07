classdef ComputeChunked < mom.ImpedanceMatrix.ComputeStrategy
    properties (Access = private)
        computeContext_ mom.ImpedanceMatrix.ComputeChunkedContext;
    end

    methods
        function obj = ComputeChunked(computeContext)
            obj.computeContext_ = computeContext;
        end
    end

    methods
        function Z = compute(obj, strategyContext, varargin)
            k = strategyContext.constants.k;
            Z0 = strategyContext.constants.Z0;
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
            rooftopPoints = strategyContext.rooftopBasisFunctions.getRooftopPoints();
            rooftopDivPoints = strategyContext.rooftopBasisFunctions.getRooftopDivPoints();

            slices = obj.computeContext_.generateSlicesFromChunks(chunks, blockSize);

            singularityMask = obj.getSingularityMask(strategyContext);
            singularityContext = mom.ImpedanceMatrix.ComputeMagPartChunked.prepareSingularityContext(strategyContext, basisFunctionCount);
            IDMask = obj.getIDMask(strategyContext);
            jacobians = obj.getJacobians(strategyContext);

            Zm = zeros(basisFunctionCount, basisFunctionCount);
            Ze = zeros(basisFunctionCount, basisFunctionCount);
            for n = 1:chunkCount
                % for debugging
                fprintf("Chunk ID=%d\n", n);

                chunk = chunks{n};
                slice = slices{n};

                rowSlice = slice.rowIndicesStart:slice.rowIndicesEnd;
                colSlice = slice.columnIndicesStart:slice.columnIndicesEnd;

                chunkContext.chunkSize = [numel(chunk.rowIndices), numel(chunk.columnIndices)];
                chunkContext.rowRooftopPoints = rooftopPoints(:, rowSlice, :);
                chunkContext.columnRooftopPoints = rooftopPoints(:, colSlice, :);
                chunkContext.singularityMaskChunk = singularityMask(chunk.rowIndices, chunk.columnIndices, :);
                chunkContext.rowRooftopDivPointsPolarity = rooftopDivPoints(:, chunk.rowIndices, p1);
                chunkContext.columnRooftopDivPointsPolarity = rooftopDivPoints(:, chunk.columnIndices, p2);

                rowQuadraturePointsPolarity = quadraturePoints(:, rowSlice, p1);
                columnQuadraturePointsPolarity = quadraturePoints(:, colSlice, p2);

                chunkContext.integralKernel = mom.ImpedanceMatrix.IntegralKernel.compute(strategyContext, rowQuadraturePointsPolarity, columnQuadraturePointsPolarity);

                [ZmBlock, singularityMaskChunkNonZeroDot] = mom.ImpedanceMatrix.ComputeMagPartChunked.computeZBlock(strategyContext, chunkContext, quadratureWeights, blockSize);
                ZeBlock = mom.ImpedanceMatrix.ComputeElPartChunked.computeZBlock(strategyContext, chunkContext, quadratureWeights, blockSize);

                Zm(chunk.rowIndices, chunk.columnIndices) = ZmBlock + mom.ImpedanceMatrix.ComputeMagPartChunked.computeMagSingularities(strategyContext, singularityContext, chunk, singularityMaskChunkNonZeroDot, IDMask);
                Ze(chunk.rowIndices, chunk.columnIndices) = ZeBlock;
            end

            Ze = Ze + mom.ImpedanceMatrix.ComputeElPartChunked.computeElSingularities(strategyContext, singularityMask, IDMask);

            Z = 1i * k * Z0 / (4 * pi) .* (Zm - (1 / k^2) .* Ze) .* jacobians;
        end
    end
end
