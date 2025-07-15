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
        function chunks = generateFull(obj, rooftops, includeDiagonalBlocks)
            [k, l] = rooftopMoM.core.impedanceMatrix.utils.generateChunkSize(obj.elementCount_, rooftops.quadrature().order());

            rooftopIndices = rooftops.indices();

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

                % Iterate over the chunk grid (p = row, q = column).
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
        end

        function chunks = generateSpecified(obj, rooftops, interactions)
        end
    end
end
