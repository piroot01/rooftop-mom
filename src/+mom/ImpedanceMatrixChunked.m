classdef ImpedanceMatrixChunked < handle
    properties
        rowChunkIndices_ uint32;
        colChunkIndices_ uint32;
    end

    methods
        function obj = ImpedanceMatrixChunked(memory_MB, order, basisFunctions)
            % input:
            %   memory_MB :: memory allocated to the matrix computation
            %   basisFunctions :: mom.BasisFunctions

            obj.generateChunkIndices(memory_MB, order, basisFunctions);
        end
    end

    methods
        function Z = compute(obj)
            chunkCount = size(obj.rowChunkIndices_, 2);

            for n = 1:chunkCount
                % for debugging
                fprintf("Chunk ID=%d, ROW=[%d, %d], COL=[%d, %d]\n", n, obj.rowChunkIndices_(1, n), obj.rowChunkIndices_(2, n), obj.colChunkIndices_(1, n), obj.colChunkIndices_(2, n));
            end
        end
    end

    methods (Access = private)
        function generateChunkIndices(obj, memory_MB, order, basisFunctions)
            % input:
            %   memory_MB :: memory allocated to the matrix computation
            %   basisFunctions :: mom.BasisFunctions

            % asssume that the distance matrix is composed out of doubles
            bytesPerElement = 8;

            function count = toElementCount(memory_MB)
                count = memory_MB * 1e6 / bytesPerElement;
            end

            function [a, b] = getChunkShape(blockCount)
                factors = [];
                for p = 1:blockCount
                    if mod(blockCount, p) == 0
                        q = blockCount / p;
                        factors = [factors; p, q];
                    end
                end

                [~, idx] = min(abs(factors(:, 1) - factors(:, 2)))
                chunkShape = factors(idx, :);

                a = chunkShape(1);
                b = chunkShape(2);
            end

            function count = getBlockCountPerChunk(elementCount, order)
                count = floor(elementCount / (order^4));
            end

            elementCount = toElementCount(memory_MB);

            if elementCount < order^4
                error("The given memory does not suffice!");
            end

            blockCountPerChunk = getBlockCountPerChunk(elementCount, order)

            [a, b] = getChunkShape(blockCountPerChunk);

            basisFuncCountPerSide = basisFunctions.getCount();

            chunkCountVert = ceil(basisFuncCountPerSide / a);
            chunkCountHorz = ceil(basisFuncCountPerSide / b);
            totalChunks = chunkCountVert * chunkCountHorz;

            rowIndices = zeros(2, totalChunks);
            colIndices = zeros(2, totalChunks);

            chunkIndex = 1;
            for m = 1:chunkCountVert
                for n = 1:chunkCountHorz
                    rStart = (m - 1) * a + 1;
                    cStart = (n - 1) * b + 1;
                    rEnd = min(m * a, basisFuncCountPerSide);
                    cEnd = min(n * b, basisFuncCountPerSide);
                    rowIndices(:, chunkIndex) = [rStart; rEnd];
                    colIndices(:, chunkIndex) = [cStart, cEnd];
                    chunkIndex = chunkIndex + 1;
                end
            end

            obj.rowChunkIndices_ = rowIndices;
            obj.colChunkIndices_ = colIndices;
        end
    end
end