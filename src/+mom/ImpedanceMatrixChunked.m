classdef ImpedanceMatrixChunked < handle
    properties (Access = private)
        quadrature_ mom.Quadrature;
        basisFunctions_ mom.BasisFunctions;

        rowChunkIndices_ uint32;
        colChunkIndices_ uint32;

        % impedance matrix constants
        minusIK0_ = -1i;
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

        function Zm = computeMag(obj)
            chunkCount = obj.getChunkCount();
            order2 = obj.quadrature_.getOrder()^2;
            weights = obj.quadrature_.getWeights();
            quadPoints = obj.basisFunctions_.getQuadPoints();
            rooftopPoints = obj.basisFunctions_.getRooftopPoints();
            basisFuncCount = obj.basisFunctions_.getCount();

            [rowStarts, rowEnds, colStarts, colEnds] = obj.generateSlices(order2);
            [r, c] = obj.generateSliceSizes();
            singularityMask = obj.generateSingularityMask();
            IDMask = obj.generateIDMask();

            Zm = zeros(basisFuncCount);
            for n = 1:chunkCount
                % for debugging
                fprintf("Chunk ID=%d, ROW=[%d, %d], COL=[%d, %d]\n", n, obj.rowChunkIndices_(1, n), obj.rowChunkIndices_(2, n), obj.colChunkIndices_(1, n), obj.colChunkIndices_(2, n));

                iSlice = rowStarts(n):rowEnds(n);
                jSlice = colStarts(n):colEnds(n);

                iQP = quadPoints(:, iSlice, :);
                jQP = quadPoints(:, jSlice, :);
                iRP = rooftopPoints(:, iSlice, :);
                jRP = rooftopPoints(:, jSlice, :);

                iChunk = obj.rowChunkIndices_(1, n):obj.rowChunkIndices_(2, n);
                jChunk = obj.colChunkIndices_(1, n):obj.colChunkIndices_(2, n);

                singularityMaskChunk = singularityMask(iChunk, jChunk, :);
                IDMaskChunk = IDMask(iChunk, jChunk);

                ZmBlock = obj.computeZmBlock(r(n), c(n), iQP, jQP, iRP, jRP, weights, order2, singularityMaskChunk, IDMaskChunk);
                Zm(iChunk, jChunk) = ZmBlock;
            end
        end

        function count = getChunkCount(obj)
            count = size(obj.rowChunkIndices_, 2);
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
                upperBound = floor(sqrt(blockCount));
                p = find(mod(blockCount, 1:upperBound) == 0);
                q = blockCount ./ p;

                [~, idx] = min(abs(p - q));
                a = p(idx);
                b = q(idx);
            end

            function count = getBlockCountPerChunk(elementCount, order)
                count = floor(elementCount / (order^4));
            end

            order = obj.quadrature_.getOrder();

            elementCount = toElementCount(memory_MB);
            if elementCount < order^4
                error('The given memory does not suffice!');
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

        function [rowStarts, rowEnds, colStarts, colEnds] = generateSlices(obj, blockSize)
            rowStarts = (obj.rowChunkIndices_(1, :) - 1) * blockSize + 1;
            rowEnds = obj.rowChunkIndices_(2, :) * blockSize;
            colStarts = (obj.colChunkIndices_(1, :) - 1) * blockSize + 1;
            colEnds = obj.colChunkIndices_(2, :) * blockSize;
        end

        function [r, c] = generateSliceSizes(obj)
            r = (obj.rowChunkIndices_(2, :) - obj.rowChunkIndices_(1, :) + 1);
            c = (obj.colChunkIndices_(2, :) - obj.colChunkIndices_(1, :) + 1);
        end

        function ZmBlock = computeZmBlock(obj, r, c, iQP, jQP, iRP, jRP, W4D, o2, SM, IDM)
            ZmBlock = zeros(r, c, 4);
            singMask = zeros(r, c, 4);
            W = repmat(W4D, r, c);

            for n = 1:obj.polarityCombCount_
                p1 = obj.polarityComb_(n, 1);
                p2 = obj.polarityComb_(n, 2);

                iQPp = iQP(:, :, p1);
                jQPp = jQP(:, :, p2);
                iRPp = iRP(:, :, p1);
                jRPp = jRP(:, :, p2);

                SMp = repelem(SM(:, :, n), o2, o2);

                D = mom.DistMatrix.getChunk(iQPp, jQPp);
                D0 = (D < eps);

                RPDot = iRPp.' * jRPp;

                RPDot0 = (RPDot < eps);
                SMpP = ~RPDot0 & SMp;
                singMask(:, :, n) = SMpP(1:o2:end, 1:o2:end);

                G = exp(-1i * D) ./ D;
                G(SMpP) = G(SMpP) - 1 ./ D(SMpP);

                I = RPDot .* G;
                I(D0) = RPDot(D0) .* obj.minusIK0_;

                I = W .* I;
                Ir = reshape(I, o2, r, o2, c);
                SS = sum(sum(Ir, 1), 3);
                ZmBlock(:, :, n) = reshape(SS, [r, c]);
            end

            ZmBlock = ZmBlock + obj.computeMagSingularities(singMask, IDM);

            ZmBlock = sum(ZmBlock, 3);
        end

        function singularities = computeMagSingularities(obj, SM, IDM)
            IDCount = obj.basisFunctions_.getIDCount();
            [RL, RLI] = obj.basisFunctions_.getRooftopLengths();
            singularities = zeros(size(SM));

            for n = 1:obj.polarityCombCount_
                singularitiesp = singularities(:, :, n);
                ss = obj.evalMagSingularities(RL, n);

                for k = 1:IDCount
                    M1 = (IDM == k);
                    M2 = M1 & SM(:, :, n);
                    [idx, ~] = find(M2);
                    ssi = RLI(:, idx);
                    singularitiesp(M2) = singularitiesp(M2) + ss(ssi).';
                end

                singularities(:, :, n) = singularitiesp;
            end
        end

        function mask = generateSingularityMask(obj)
            basisFuncsCount = obj.basisFunctions_.getCount();
            basisFuncs = obj.basisFunctions_.getFaces();

            mask = false(basisFuncsCount);

            for n = 1:obj.polarityCombCount_
                p1 = obj.polarityComb_(n, 1);
                p2 = obj.polarityComb_(n, 2);

                mask(:, :, n) = (basisFuncs(p1, :).' == basisFuncs(p2, :));
            end
        end

        function mask = generateIDMask(obj)
            basisFuncIDs = double(obj.basisFunctions_.getFuncIDs());
            mask = (basisFuncIDs == basisFuncIDs.') .* (basisFuncIDs.' * ones(1, numel(basisFuncIDs)));
        end

        function singularities = evalMagSingularities(~, xs, polarity)
            if (polarity == 1 || polarity == 4)
                singularities = mom.singularity.singularity_1(xs(1, :), xs(2, :));
            else
                singularities = mom.singularity.singularity_5(xs(1, :), xs(2, :));
            end
        end
    end
end