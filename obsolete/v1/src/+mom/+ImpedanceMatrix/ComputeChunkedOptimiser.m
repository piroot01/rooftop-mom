classdef ComputeChunkedOptimiser < handle
    properties
        computeStrategy_;
        computeChunkedContext_ mom.ImpedanceMatrix.ComputeChunkedContext;
    end

    methods
        function obj = ComputeChunkedOptimiser(computeStrategy, computeChunkedContext)
            obj.computeStrategy_ = computeStrategy;
            obj.computeChunkedContext_ = computeChunkedContext;
        end
    end

    methods
        % NOTE: does not work with multiple antennas yet
        function Z = compute(obj, strategyContext, gridSizes)
            computeStrategyOptimiser = mom.ImpedanceMatrix.ComputeStrategyOptimiser(strategyContext.basisFunctions, strategyContext.rooftopBasisFunctions.getUniqueUniqueRooftopSizes(strategyContext.basisFunctions.getIDs()), gridSizes);
            [masksCell, interactionsCell] = computeStrategyOptimiser.getAllMasksAndInteractions();

            [mask, interactions, chunks] = obj.generateMaskAndInteractionsAndChunks(strategyContext, masksCell, interactionsCell);

            ZUniqueValues = obj.computeStrategy_.compute(strategyContext, chunks);

            Z = computeStrategyOptimiser.applyMask(ZUniqueValues, mask, interactions);
        end
    end

    methods (Access = private)
        function [mask, interactions, chunks] = generateMaskAndInteractionsAndChunks(obj, strategyContext, masksCell, interactionsCell)
            basisFunctions = strategyContext.basisFunctions.getFacesAsCell();
            IDCount = strategyContext.basisFunctions.getIDCount();
            IDMask = obj.computeStrategy_.getIDMask(strategyContext);
            quadratureOrder = strategyContext.quadrature.getOrder();

            basisFunctionCounts = cellfun(@(n) size(n, 2), basisFunctions);

            maskMaxima = cellfun(@(m) max(m(:)), masksCell);
            maskOffsets = [0, cumsum(maskMaxima(1:(end - 1)))];

            chunkOffsets = [0, cumsum(basisFunctionCounts(1:(end - 1)))];

            interactionCountPerCell = cellfun(@numel, interactionsCell);
            interactionOffsets = [0, cumsum(interactionCountPerCell(1:(end - 1)))];

            mask = zeros(size(IDMask), 'uint32');
            interactions = zeros(sum(interactionCountPerCell), 1);
            chunksDiagonalSubmatrices = cell(1, IDCount);

            for n = 1:IDCount
                chunks = obj.computeChunkedContext_.generateChunksSpecificInteractions(quadratureOrder, interactionsCell{n}, basisFunctionCounts(n));

                chunkOffset = chunkOffsets(n);
                for k = 1:numel(chunks)
                    chunks{k}.rowIndices = chunks{k}.rowIndices + chunkOffset;
                    chunks{k}.columnIndices = chunks{k}.columnIndices + chunkOffset;
                end
                chunksDiagonalSubmatrices{n} = chunks;

                translatedInteractions = obj.translateLinearIndices(size(mask), [basisFunctionCounts(n), basisFunctionCounts(n)], interactionsCell{n}, chunkOffset, chunkOffset);

                startIndex = interactionOffsets(n) + 1;
                endIndex = interactionOffsets(n) + interactionCountPerCell(n);
                interactions(startIndex:endIndex) = translatedInteractions;

                mask(IDMask == n) = masksCell{n} + maskOffsets(n);
            end

            if IDCount > 1
                chunksOffdiagonalSubmatrices = obj.computeChunkedContext_.generateChunksAllInteractions(quadratureOrder, basisFunctionCounts);
            else
                chunksOffdiagonalSubmatrices = {};
            end

            chunks = [[chunksDiagonalSubmatrices{:}], chunksOffdiagonalSubmatrices];
        end

        function indices = translateLinearIndices(~, maskSize, blockSize, blockIndices, rowStart, colStart)
            [rows, cols] = ind2sub(blockSize, blockIndices);
            rows = rows + rowStart;
            cols = cols + colStart;
            indices = uint32(sub2ind(maskSize, rows, cols));
        end
    end
end
