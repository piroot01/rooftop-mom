classdef ComputeMagPartChunked < mom.ImpedanceMatrix.ComputeStrategy
    properties (Access = private)
        computeContext_ mom.ImpedanceMatrix.ComputeChunkedContext;
    end

    methods
        function obj = ComputeMagPartChunked(computeContext)
            obj.computeContext_ = computeContext;
        end
    end

    methods
        function Z = compute(obj, strategyContext)
            chunkCount = obj.computeContext_.getChunkCount();
            order2 = strategyContext.quadrature.getOrder()^2;
            weights = strategyContext.quadrature.getWeights();
            basisFuncCount = strategyContext.basisFunctions.getCount();

            quadPoints = strategyContext.basisFunctions.getQuadPoints();
            rooftopPoints = strategyContext.basisFunctions.getRooftopPoints();

            [rowStarts, rowEnds, colStarts, colEnds] = obj.computeContext_.getSlices(order2);
            [r, c] = obj.computeContext_.getSliceSizes();

            singularityMask = obj.getSingularityMask(strategyContext);
            IDMask = obj.getIDMask(strategyContext);
            jacobians = obj.getJacobians(strategyContext);

            Z = zeros(basisFuncCount);
            for n = 1:chunkCount
                % for debugging
                fprintf("Chunk ID=%d\n", n);

                iSlice = rowStarts(n):rowEnds(n);
                jSlice = colStarts(n):colEnds(n);

                iQP = quadPoints(:, iSlice, :);
                jQP = quadPoints(:, jSlice, :);
                iRP = rooftopPoints(:, iSlice, :);
                jRP = rooftopPoints(:, jSlice, :);

                iChunk = obj.computeContext_.rowChunkIndices(1, n):obj.computeContext_.rowChunkIndices(2, n);
                jChunk = obj.computeContext_.colChunkIndices(1, n):obj.computeContext_.colChunkIndices(2, n);

                singularityMaskChunk = singularityMask(iChunk, jChunk, :);
                IDMaskChunk = IDMask(iChunk, jChunk);

                ZBlock = obj.computeZBlock(strategyContext, r(n), c(n), iQP, jQP, iRP, jRP, weights, order2, singularityMaskChunk, IDMaskChunk);
                Z(iChunk, jChunk) = ZBlock;
            end

            Z = Z .* jacobians;
        end
    end

    methods (Access = private)
        function ZBlock = computeZBlock(obj, strategyContext, r, c, iQP, jQP, iRP, jRP, W4D, o2, SM, IDM)
            ZBlock4 = zeros(r, c, 4);
            singMask = false(r, c, 4);
            W = repmat(W4D, r, c);

            for n = 1:strategyContext.constants.polaritiesCount
                p1 = strategyContext.constants.polarities(n, 1);
                p2 = strategyContext.constants.polarities(n, 2);

                iQPp = iQP(:, :, p1);
                jQPp = jQP(:, :, p2);
                iRPp = iRP(:, :, p1);
                jRPp = jRP(:, :, p2);

                SMp = repelem(SM(:, :, n), o2, o2);

                D = mom.DistMatrix.getChunk(iQPp, jQPp);
                D0 = (D < eps);

                RPDot = iRPp.' * jRPp;
                RPDot0 = (RPDot < eps);

                SMpRP = ~RPDot0 & SMp;
                singMask(:, :, n) = SMpRP(1:o2:end, 1:o2:end);

                G = exp(-1i * D) ./ D;
                G(SMpRP) = G(SMpRP) - 1 ./ D(SMpRP);

                I = RPDot .* G;
                I(D0) = RPDot(D0) .* -1i;

                IW = I .* W;
                Ir = reshape(IW, o2, r, o2, c);
                SS = sum(sum(Ir, 1), 3);
                ZBlock4(:, :, n) = reshape(SS, r, c);
            end

            ZBlock4 = ZBlock4 + obj.computeMagSingularities(strategyContext, singMask, IDM);
            ZBlock = sum(ZBlock4, 3);
        end

        function singularities = computeMagSingularities(obj, strategyContext, SM, IDM)
            IDCount = strategyContext.basisFunctions.getIDCount();
            [RL, RLI] = strategyContext.basisFunctions.getRooftopLengths();
            singularities = zeros(size(SM));

            for n = 1:strategyContext.constants.polaritiesCount
                SP = singularities(:, :, n);
                ss = obj.evalMagSingularities(RL, n);

                for k = 1:IDCount
                    M1 = (IDM == k);
                    M2 = M1 & SM(:, :, n);
                    [idx, ~] = find(M2);
                    ssi = RLI(:, idx);
                    SP(M2) = SP(M2) + ss(ssi).';
                end

                singularities(:, :, n) = SP;
            end
        end

        function singularities = evalMagSingularities(~, xs, polarity)
            if (polarity == 1 || polarity == 4)
                singularities = mom.ImpedanceMatrix.singularity.type1(xs(1, :), xs(2, :));
            else
                singularities = mom.ImpedanceMatrix.singularity.type2(xs(1, :), xs(2, :));
            end
        end
    end
end