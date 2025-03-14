classdef ComputeChunkedContext < handle
    properties (Access = public)
        rowChunkIndices uint32;
        colChunkIndices uint32;
    end

    methods
        function obj = ComputeChunkedContext(memory_MB, quadOrder, basisFuncCount)
            obj.generateChunkIndices(memory_MB, quadOrder, basisFuncCount);
        end
    end

    methods
        function count = getChunkCount(obj)
            count = size(obj.rowChunkIndices, 2);
        end

        function [rowStarts, rowEnds, colStarts, colEnds] = getSlices(obj, blockSize)
            rowStarts = (obj.rowChunkIndices(1, :) - 1) * blockSize + 1;
            rowEnds = obj.rowChunkIndices(2, :) * blockSize;
            colStarts = (obj.colChunkIndices(1, :) - 1) * blockSize + 1;
            colEnds = obj.colChunkIndices(2, :) * blockSize;
        end

        function [r, c] = getSliceSizes(obj)
            r = (obj.rowChunkIndices(2, :) - obj.rowChunkIndices(1, :) + 1);
            c = (obj.colChunkIndices(2, :) - obj.colChunkIndices(1, :) + 1);
        end
    end

    methods (Access = private)
        function generateChunkIndices(obj, memory_MB, quadOrder, basisFuncCount)
            bytesPerElement = 8;
            bytesInMegabyte = 1e6;

            function count = toElementCount(memory_MB)
                count = memory_MB * bytesInMegabyte / bytesPerElement;
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