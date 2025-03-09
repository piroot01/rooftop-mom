classdef ImpedanceMatrixChunked < handle
    properties (Access = private)
        quadrature_ mom.Quadrature;
        basisFunctions_ mom.BasisFunctions;

        rowChunkIndices_ uint32;
        colChunkIndices_ uint32;
    end

    properties (Access = private, Constant)
        polarityComb_ = [0, 0; 0, 1; 1, 0; 1, 1] + 1;
        polarityCombCount_ = 4;
    end

    methods
        function obj = ImpedanceMatrixChunked(quadrature, basisFunctions, memory_MB)
            % input:
            %   quadrature :: mom.Quadrature
            %   basisFunctions :: mom.BasisFunctions
            %   memory_MB :: memory allocated to the matrix computation
            obj.quadrature_ = quadrature;
            obj.basisFunctions_ = basisFunctions;

            obj.generateChunkIndices(memory_MB);
        end
    end

    methods
        function Z = compute(obj)
            chunkCount = size(obj.rowChunkIndices_, 2);

            quadPoints = obj.basisFunctions_.getQuadPoints();
            rooftopPoints = obj.basisFunctions_.getRooftopPoints();

            weights4D = obj.quadrature_.getWeights();
            order2 = obj.quadrature_.getOrder()^2;

            for n = 1:chunkCount
                % for debugging
                fprintf("Chunk ID=%d, ROW=[%d, %d], COL=[%d, %d]\n", n, obj.rowChunkIndices_(1, n), obj.rowChunkIndices_(2, n), obj.colChunkIndices_(1, n), obj.colChunkIndices_(2, n));

                iStart = obj.rowChunkIndices_(1, n);
                iEnd = obj.rowChunkIndices_(2, n);
                jStart = obj.colChunkIndices_(1, n);
                jEnd = obj.colChunkIndices_(2, n);

                iSlice = ((iStart - 1) * order2 + 1):(iEnd * order2);
                jSlice = ((jStart - 1) * order2 + 1):(jEnd * order2);

                r = numel(iSlice) / order2;
                s = numel(jSlice) / order2;

                iQP = quadPoints(:, iSlice, :);
                jQP = quadPoints(:, jSlice, :);

                iRP = rooftopPoints(:, iSlice, :);
                jRP = rooftopPoints(:, jSlice, :);

                ZBlock = obj.computeBlock(r, s, iQP, jQP, iRP, jRP, weights4D, order2);
            end
        end
    end

    methods (Access = private)
        function generateChunkIndices(obj, memory_MB)
            bytesPerElement = 8;
            bytesInMegabyte = 1e6;

            function count = toElementCount(memory_MB)
                count = memory_MB * bytesInMegabyte / bytesPerElement;
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

            order = obj.quadrature_.getOrder();

            elementCount = toElementCount(memory_MB);
            if elementCount < order^4
                error("The given memory does not suffice!");
            end

            blockCountPerChunk = getBlockCountPerChunk(elementCount, order);
            [a, b] = getChunkShape(blockCountPerChunk);
            basisFuncCountPerSide = obj.basisFunctions_.getCount();

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

        function Zblock = computeBlock(obj, m, n, iQP, jQP, iRP, jRP, W4D, o2)
            Zblock = zeros(m, n);
            W = repmat(W4D, m, n);
            for n = 1
                iQPp = iQP(:, :, obj.polarityComb_(n, 1))
                jQPp = jQP(:, :, obj.polarityComb_(n, 2))
                iRPp = iRP(:, :, obj.polarityComb_(n, 1))
                jRPp = jRP(:, :, obj.polarityComb_(n, 2))

                D = mom.DistMatrix.getChunk(iQPp, jQPp)
                RPDot = iRPp.' * jRPp
                G = exp(-1i * D) ./ D;
                G(D == 0) = 0;
                I = W .* RPDot .* G
                Ir = reshape(I, o2, m, o2, []);
                Zblock = Zblock + squeeze(sum(sum(Ir, 1), 3));
            end

            disp(Zblock)
        end
    end
end