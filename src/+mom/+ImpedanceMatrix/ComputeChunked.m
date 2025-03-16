classdef ComputeChunked < mom.ImpedanceMatrix.ComputeStrategy
    properties (Access = private)
        computeContext_ mom.ImpedanceMatrix.ComputeChunkedContext;
    end

    methods
        function obj = ComputeChunked(computeContext)
            obj.computeContext_ = computeContext;
        end
    end

    methods
        function Z = compute(obj, strategyContext)
            k = strategyContext.constants.k;
            Z0 = strategyContext.constants.Z0;

            chunkCount = obj.computeContext_.getChunkCount();
            order2 = strategyContext.quadrature.getOrder()^2;
            weights = strategyContext.quadrature.getWeights();
            basisFuncCount = strategyContext.basisFunctions.getCount();

            quadPoints = strategyContext.basisFunctions.getQuadPoints();
            rooftopPoints = strategyContext.basisFunctions.getRooftopPoints();
            rooftopPointsDiv = strategyContext.basisFunctions.getRooftopPointsDiv();

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
                iChunk = obj.computeContext_.rowChunkIndices(1, n):obj.computeContext_.rowChunkIndices(2, n);
                jChunk = obj.computeContext_.colChunkIndices(1, n):obj.computeContext_.colChunkIndices(2, n);

                iQP = quadPoints(:, iSlice, :);
                jQP = quadPoints(:, jSlice, :);
                iRP = rooftopPoints(:, iSlice, :);
                jRP = rooftopPoints(:, jSlice, :);
                iRPDiv = rooftopPointsDiv(:, iChunk, :);
                jRPDiv = rooftopPointsDiv(:, jChunk, :);

                singularityMaskChunk = singularityMask(iChunk, jChunk, :);
                IDMaskChunk = IDMask(iChunk, jChunk);

                ZBlock = obj.computeZBlock(strategyContext, r(n), c(n), ...
                    iQP, jQP, iRP, jRP, iRPDiv, jRPDiv, ...
                    weights, order2, singularityMaskChunk, IDMaskChunk);

                Z(iChunk, jChunk) = ZBlock;
            end

            Z = 1i * k * Z0 / (4 * pi) .* Z .* jacobians;
        end
    end

    methods (Access = private)
        function ZBlock = computeZBlock(obj, strategyContext, r, c, iQP, jQP, iRP, jRP, iRPDiv, jRPDiv, W4D, o2, SM, IDM)
            k = strategyContext.constants.k;

            ZmBlock4 = zeros(r, c, 4);
            ZeBlock4 = zeros(r, c, 4);
            RPDiv = zeros(r, c, 4);
            singMaskMag = false(r, c, 4);
            W = repmat(W4D, r, c);

            for n = 1:strategyContext.constants.polaritiesCount
                p1 = strategyContext.constants.polarities(n, 1);
                p2 = strategyContext.constants.polarities(n, 2);

                iQPp = iQP(:, :, p1);
                jQPp = jQP(:, :, p2);

                D = mom.DistMatrix.getChunk(iQPp, jQPp);
                D0 = (D < eps);

                G = exp(-1i * k * D) ./ D;

                % magnetical part
                iRPp = iRP(:, :, p1);
                jRPp = jRP(:, :, p2);
                RPDot = iRPp.' * jRPp;
                RPDot0 = (RPDot < eps);

                SMp = repelem(SM(:, :, n), o2, o2);
                SMpRP = ~RPDot0 & SMp;
                singMaskMag(:, :, n) = SMpRP(1:o2:end, 1:o2:end);

                Gm = G;
                Gm(SMpRP) = Gm(SMpRP) - 1 ./ D(SMpRP);

                Im = RPDot .* Gm;
                Im(D0) = RPDot(D0) .* -1i * k;
                ImW = Im .* W;
                Imr = reshape(ImW, o2, r, o2, c);
                SSm = sum(sum(Imr, 1), 3);
                ZmBlock4(:, :, n) = reshape(SSm, r, c);

                % electrical part
                iRPDivp = iRPDiv(:, :, p1);
                jRPDivp = jRPDiv(:, :, p2);
                RPDiv(:, :, n) = iRPDivp.' * jRPDivp;

                Ge = G;
                Ge(SMp) = Ge(SMp) - 1 ./ D(SMp);

                Ie = Ge;
                Ie(D0) = -1i * k;
                IeW = Ie .* W;
                Ier = reshape(IeW, o2, r, o2, c);
                SSe = sum(sum(Ier, 1), 3);
                ZeBlock4(:, :, n) = reshape(SSe, r, c);
            end

            ZmSings4 = obj.computeMagSingularities(strategyContext, singMaskMag, IDM);
            ZeSings4 = obj.computeElSingularities(strategyContext, SM, IDM);

            ZBlock4 = (ZmBlock4 + ZmSings4) - 1 / k^2 * ((ZeBlock4 + ZeSings4) .* RPDiv);
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

        function singularities = computeElSingularities(obj, strategyContext, SM, IDM)
            IDCount = strategyContext.basisFunctions.getIDCount();
            [RL, RLI] = strategyContext.basisFunctions.getRooftopLengths();
            singularities = zeros(size(SM));

            ss = obj.evalElSingularities(RL);

            for n = 1:strategyContext.constants.polaritiesCount
                SP = singularities(:, :, n);
                M2 = SM(:, :, n);

                for k = 1:IDCount
                    M1 = (IDM == k);
                    combinedMask = M1 & M2;

                    if any(combinedMask(:))
                        [idx, ~] = find(M2);
                        ssi = RLI(:, idx);
                        SP(combinedMask) = SP(combinedMask) + ss(ssi).';
                    end
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
        
        function singularities = evalElSingularities(~, xs)
            singularities = mom.ImpedanceMatrix.singularity.type3(xs(1, :), xs(2, :));
        end
    end
end