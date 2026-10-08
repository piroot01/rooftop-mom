classdef AssemblerChunkPlanner < handle
    properties (Access = private)
        elementCount_ uint32;
    end

    methods
        function obj = AssemblerChunkPlanner(elementCount)
            obj.elementCount_ = elementCount;
        end
    end

    methods
        function chunks = generateFull(obj, context, includeDiagonalBlocks)
            TRACE("AssemblerChunkPlanner.generateFull: Initialized the chunk generation procedure.");

            [k, l] = rooftopMoM.core.impedanceMatrix.utils.generateChunkSize(obj.elementCount_, context.quadrature().order());

            rooftopIndices = context.rooftops().indices();

            [offsets, sizes] = rooftopMoM.core.impedanceMatrix.utils.generateSubmatrixOffsetsAndSizes(rooftopIndices);
            indices = rooftopMoM.core.impedanceMatrix.utils.generateSubmatrixIndices(rooftopIndices, includeDiagonalBlocks);

            % Vectorize the calculation of how many chunks are needed for all specified submatrices.
            numChunksInRows = ceil(sizes(1, indices) ./ k);
            numChunksInCols = ceil(sizes(2, indices) ./ l);
            totalChunks = sum(numChunksInRows .* numChunksInCols);

            % Pre-allocate the cell array to prevent resizing in the loop.
            chunks = cell(1, totalChunks);
            chunkIndex = 1;

            % Iterate over each specified submatrix index.
            for m = indices
                % Size of the current submatrix.
                sz = sizes(:, m);

                % Offset of the submatrix.
                off = offsets(:, m);

                % Determine the grid dimensions for the chunks within this submatrix.
                rowCount = numChunksInRows(m);
                columnCount = numChunksInCols(m);

                % Iterate over the chunk grid.
                for p = 1:rowCount
                    for q = 1:columnCount
                        % Calculate chunk boundaries, local to the submatrix.
                        rowStart = (p - 1) * k + 1;
                        rowEnd   = min(p * k, sz(1));
                        colStart = (q - 1) * l + 1;
                        colEnd   = min(q * l, sz(2));

                        % Create a struct for the current chunk.
                        chunk.rowIndices    = (rowStart:rowEnd) + off(1);
                        chunk.columnIndices = (colStart:colEnd) + off(2);

                        chunks{chunkIndex} = chunk;
                        chunkIndex = chunkIndex + 1;
                    end
                end
            end

            TRACE("AssemblerChunkPlanner.generateFull: Finished chunk generation, total of %d chunks.", totalChunks);
        end

        function chunks = generateSpecified(obj, context, interactions, count)
            TRACE("AssemblerChunkPlanner.generateSpecified: Initialized the chunk generation procedure.");

            blockCount = rooftopMoM.core.impedanceMatrix.utils.generateBlockCountPerChunk(obj.elementCount_, context.quadrature().order());
            groups = rooftopMoM.core.impedanceMatrix.utils.generateInteractionGroupsByColumn(interactions, count);

            % Calculate the number of chunks needed for each interaction group.
            numChunksInGroups = cellfun(@(group) ceil(numel(group) / blockCount), groups);
            totalChunks = sum(numChunksInGroups);

            % Pre-allocate the output cell array to its final size to avoid slow resizing.
            chunks = cell(1, totalChunks);
            chunkIndex = 1;

            % Iterate through each group of interactions.
            for m = 1:numel(groups)
                currentGroup = groups{m};
                groupSize = numel(currentGroup);

                % Iterate through the chunks required for the current group.
                for n = 1:numChunksInGroups(m)
                    % Determine the slice of indices from the group for this specific chunk.
                    startIdx = (n - 1) * blockCount + 1;
                    endIdx = min(n * blockCount, groupSize);

                    % Get the actual linear indices of the matrix elements for this chunk.
                    chunkIndices = currentGroup(startIdx:endIdx);

                    % Convert the linear indices to subscripts.
                    % The full matrix is assumed to be square with size [count x count].
                    [rows, cols] = ind2sub([count, count], chunkIndices);

                    % Define a rectangular chunk.
                    chunk.rowIndices = rows(1):rows(end);
                    chunk.columnIndices = cols(1):cols(end);

                    chunks{chunkIndex} = chunk;
                    chunkIndex = chunkIndex + 1;
                end
            end

            TRACE("AssemblerChunkPlanner.generateFull: Finished chunk generation, total of %d chunks.", totalChunks);
        end
    end
end
