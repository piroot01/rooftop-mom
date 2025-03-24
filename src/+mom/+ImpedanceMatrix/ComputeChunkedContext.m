classdef ComputeChunkedContext < handle
    properties (Access = public)
        memory_MB_ double;

        rowChunkIndices uint32;
        colChunkIndices uint32;
    end

    properties (Access = private, Constant)
        bytesPerElement_ = double(8);
        bytesInMegabyte_ = double(1e6);
    end

    methods
        function obj = ComputeChunkedContext(memory_MB)
            obj.memory_MB_ = double(memory_MB);

            fprintf("Memory used: %lu MB\n", obj.memory_MB_);
        end
    end

    methods
        function chunks = generateChunksAllInteractions(obj, quadratureOrder, basisFunctionCounts)
            elementCount = obj.memory_MB_ * obj.bytesInMegabyte_ / obj.bytesPerElement_;
            if elementCount < quadratureOrder^4
                error('The given memory does not suffice!');
            end

            function [a, b] = getChunkShape(blockCountPerChunk)
                upperBound = floor(sqrt(blockCountPerChunk));
                p = find(mod(blockCountPerChunk, 1:upperBound) == 0);
                q = blockCountPerChunk ./ p;

                [~, idx] = min(abs(p - q));
                a = p(idx);
                b = q(idx);
            end

            elementCountPerBlock = quadratureOrder^4;
            blockCountPerChunk = floor(double(elementCount) / elementCountPerBlock);
            [a, b] = getChunkShape(blockCountPerChunk);

            diagonalSubmatrixCount = numel(basisFunctionCounts);
            diagonalSubmatrixOffsets = [0, cumsum(basisFunctionCounts(1:(end - 1)))];

            [rowIndices, colIndices] = ndgrid(1:diagonalSubmatrixCount);
            submatrixOffsets = [diagonalSubmatrixOffsets(rowIndices(:)); diagonalSubmatrixOffsets(colIndices(:))];
            submatrixSizes = [basisFunctionCounts(rowIndices(:)); basisFunctionCounts(colIndices(:))];

            if diagonalSubmatrixCount > 1
                allSubmatrixIndices = 1:size(submatrixSizes, 2);
                allSubmatrixIndicesCount = numel(allSubmatrixIndices);
                skipDiagonalSubmatrixIndices = 1:(diagonalSubmatrixCount + 1):allSubmatrixIndicesCount;
                submatrixIndices = setdiff(allSubmatrixIndices, skipDiagonalSubmatrixIndices);
            else
                submatrixIndices = 1;
            end

            chunks = cell(1);
            chunkIndex = 1;

            for m = submatrixIndices
                submatrixSize = submatrixSizes(:, m);
                submatrixOffset = submatrixOffsets(:, m);

                rowChunkCount = ceil(submatrixSize(1) / a);
                columnChunkCount = ceil(submatrixSize(2) / b);

                [rowIndices, colIndices] = meshgrid(1:rowChunkCount, 1:columnChunkCount);
                indices = [rowIndices(:), colIndices(:)].';

                for n = indices
                    p = n(1);
                    q = n(2);

                    rowStart = (p - 1) * a + 1;
                    colStart = (q - 1) * b + 1;
                    rowEnd = min(p * a, submatrixSize(1));
                    colEnd = min(q * b, submatrixSize(2));

                    chunk.rowIndices = (rowStart:rowEnd) + submatrixOffset(1);
                    chunk.columnIndices = (colStart:colEnd) + submatrixOffset(2);
                    chunks{chunkIndex} = chunk;
                    chunkIndex = chunkIndex + 1;
                end
            end
        end

        function chunks = generateChunksSpecificInteractions(obj, quadratureOrder, interactions, basisFunctionCount)
            elementCount = obj.memory_MB_ * obj.bytesInMegabyte_ / obj.bytesPerElement_;
            if elementCount < quadratureOrder^4
                error('The given memory does not suffice!');
            end

            elementCountPerBlock = quadratureOrder^4;
            blockCountPerChunk = floor(elementCount / elementCountPerBlock);
            % interactionCount = numel(interactions);

            sortedInteractions = sort(interactions).';

            columnIndices = floor(double((sortedInteractions - 1)) / basisFunctionCount) + 1;
            splitCols = [diff(columnIndices) ~= 0, true];
            splitColIndices = find(splitCols);
            groupSizes = diff([0, splitColIndices]);
            chunkGroups = mat2cell(sortedInteractions(:), groupSizes);

            %{
            groupBoundaries = [true, diff(sortedInteractions) ~= 1];
            groupIndices = cumsum(groupBoundaries)
            groupStarts = find(groupBoundaries);
            groupMask = (1:interactionCount) - groupStarts(groupIndices) + 1
            groupMaskFiltered = (mod(groupMask - 1, basisFunctionCount + 1) == 0) & (groupMask > 1)
            splitGroupMask = groupBoundaries | groupMaskFiltered;
            splitGroupMaskIndices = find(splitGroupMask);
            groupSizes = diff([splitGroupMaskIndices, interactionCount + 1])
            chunkGrops = mat2cell(sortedInteractions(:), groupSizes)
            %}

            groupCount = numel(groupSizes);

            chunks = cell(1);
            chunkIndex = 1;
            for m = 1:groupCount
                chunkCount = ceil(numel(chunkGroups{m}) / double(blockCountPerChunk));
                for n = 1:chunkCount
                    slice = ((n - 1) * blockCountPerChunk + 1):min(numel(chunkGroups{m}),(n * blockCountPerChunk));
                    chunkIndices = chunkGroups{m}(slice);
                    [rows, cols] = ind2sub([basisFunctionCount, basisFunctionCount], chunkIndices);
                    rowStart = rows(1);
                    colStart = cols(1);
                    rowEnd = rows(end);
                    colEnd = cols(end);

                    chunk.rowIndices = rowStart:rowEnd;
                    chunk.columnIndices = colStart:colEnd;
                    chunks{chunkIndex} = chunk;
                    chunkIndex = chunkIndex + 1;
                end
            end
        end

        function slices = generateSlicesFromChunks(~, chunks, blockSize)
            chunkCount = numel(chunks);
            slices = cell(chunkCount, 1);

            for m = 1:chunkCount
                chunk = chunks{m};

                rowMin = chunk.rowIndices(1);
                rowMax = chunk.rowIndices(end);
                colMin = chunk.columnIndices(1);
                colMax = chunk.columnIndices(end);

                rowStart = (rowMin - 1) * blockSize + 1;
                rowEnd = rowMax * blockSize;
                colStart = (colMin - 1) * blockSize + 1;
                colEnd = colMax * blockSize;

                sliceStruct = struct(...
                    'rowIndicesStart', rowStart, ...
                    'rowIndicesEnd', rowEnd, ...
                    'columnIndicesStart', colStart, ...
                    'columnIndicesEnd', colEnd ...
                );

                slices{m} = sliceStruct;
            end
        end
    end

    methods (Access = private)
        function generateChunkIndices(obj, memory_MB, quadOrder, basisFuncCount)
            function count = toElementCount(obj, memory_MB)
                count = memory_MB * obj.bytesInMegabyte_ / obj.bytesPerElement_;
            end

            function [a, b] = getChunkShape(blockCount)
                upperBound = floor(sqrt(blockCount));
                p = find(mod(blockCount, 1:upperBound) == 0);
                q = blockCount ./ p;

                [~, idx] = min(abs(p - q));
                a = p(idx);
                b = q(idx);
            end

            function count = getBlockCountPerChunk(elementCount, quadOrder)
                count = floor(elementCount / (quadOrder^4));
            end

            elementCount = toElementCount(memory_MB);
            if elementCount < quadOrder^4
                error('The given memory does not suffice!');
            end

            blockCountPerChunk = getBlockCountPerChunk(elementCount, quadOrder);
            [a, b] = getChunkShape(blockCountPerChunk);

            chunkCountVert = ceil(basisFuncCount / a);
            chunkCountHorz = ceil(basisFuncCount / b);
            totalChunks = chunkCountVert * chunkCountHorz;

            rowIndices = zeros(2, totalChunks);
            colIndices = zeros(2, totalChunks);

            chunkIndex = 1;
            for m = 1:chunkCountVert
                for n = 1:chunkCountHorz
                    rStart = (m - 1) * a + 1;
                    cStart = (n - 1) * b + 1;
                    rEnd = min(m * a, basisFuncCount);
                    cEnd = min(n * b, basisFuncCount);
                    rowIndices(:, chunkIndex) = [rStart; rEnd];
                    colIndices(:, chunkIndex) = [cStart, cEnd];
                    chunkIndex = chunkIndex + 1;
                end
            end

            obj.rowChunkIndices = rowIndices;
            obj.colChunkIndices = colIndices;
        end
    end
end
