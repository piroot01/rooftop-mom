classdef InteractionSymmetriesAccelerator < rooftopMoM.core.impedanceMatrix.InteractionSymmetriesAcceleratorImpl
    properties (Access = private)
        strategy_;
    end

    methods
        function obj = InteractionSymmetriesAccelerator(strategy)
            obj.strategy_ = strategy;
        end
    end

    methods
        function Z = generateFull(obj, context)
            [maskCell, interactionIndicesCell] = obj.generateInteractions(context);

            maskCell{1}
            interactionIndicesCell{1}

            [mask, interactionIndices, chunks] = obj.generateGlobal_(context, maskCell, interactionIndicesCell);

            values = obj.strategy_.generateSpecified(context, chunks);

            Z = rooftopMoM.core.impedanceMatrix.utils.fillMatrixUsingMaskAndInteractionIndices(values, mask, interactionIndices);
        end
    end

    methods (Access = private)
        function [mask, interactionIndices, chunks] = generateGlobal_(obj, context, maskCell, interactionIndicesCell)
            chunkPlanner = obj.strategy_.chunkPlanner();

            % Copute the mask offsets.
            maskMaxima = cellfun(@(m) max(m(:)), maskCell);
            maskOffsets = [0, cumsum(maskMaxima(1:(end - 1)))];

            % Compute the interaction indices offsets
            interactionIndicesCounts = cellfun(@numel, interactionIndicesCell);
            interactionIndicesOffsets = [0, cumsum(interactionIndicesCounts(1:(end - 1)))];

            % Compute the chunk offsets
            basisSizes = context.rooftops().counts();
            chunkOffsets = [0, cumsum(basisSizes(1:(end - 1)))];

            basisSize = context.rooftops().count();
            meshCount = context.rooftops().mesh().count();

            mask = zeros(basisSize, basisSize, 'uint32');
            interactionIndices = zeros(sum(interactionIndicesCounts), 1, 'uint32');
            chunksDiagonalSubmatrices = cell(1, meshCount);

            % TODO: Extract to standalone function
            applyChunkOffset = @(chunk, offset) struct( ...
                'rowIndices', chunk.rowIndices + offset, ...
                'columnIndices', chunk.columnIndices + offset ...
            );

            for n = 1:meshCount
                % Shift and inset the mask.
                startIdx = chunkOffsets(n) + 1;
                endIdx = chunkOffsets(n) + basisSizes(n);
                range = startIdx:endIdx;

                mask(range, range) = rooftopMoM.core.impedanceMatrix.utils.shiftIndices(maskCell{n}, maskOffsets(n));

                % Offset the interaction indices.
                localInteractionIndices = rooftopMoM.core.impedanceMatrix.utils.transformIndices( ...
                    interactionIndicesCell{n}, [basisSize, basisSize], size(maskCell{n}), chunkOffsets(n), chunkOffsets(n));

                % Insert the interaction
                startIdx = interactionIndicesOffsets(n) + 1;
                endIdx = interactionIndicesOffsets(n) + interactionIndicesCounts(n);
                interactionIndices(startIdx:endIdx) = localInteractionIndices;

                % Shift the chunks.
                localChunks = chunkPlanner.generateSpecified(context, interactionIndicesCell{n}, basisSizes(n));

                % Shift all chunk indices by the current offset in one go.
                localChunks = cellfun(@(chunk) applyChunkOffset(chunk, chunkOffsets(n)), localChunks, 'UniformOutput', false);

                chunksDiagonalSubmatrices{n} = localChunks;
            end

            % TODO: Determine all chunks beside the diagonal ones.
            chunks = chunksDiagonalSubmatrices{:};
        end
    end
end
