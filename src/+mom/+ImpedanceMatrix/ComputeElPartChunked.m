classdef ComputeElPartChunked < mom.ImpedanceMatrix.ComputeStrategy
    properties (Access = private)
        computeContext_ mom.ImpedanceMatrix.ComputeChunkedContext;
    end

    methods
        function obj = ComputeElPartChunked(computeContext)
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
                iRPDiv = rooftopPointsDiv(:, iChunk, :);
                jRPDiv = rooftopPointsDiv(:, jChunk, :);

                singularityMaskChunk = singularityMask(iChunk, jChunk, :);
                IDMaskChunk = IDMask(iChunk, jChunk);

                ZBlock = obj.computeZBlock(strategyContext, r(n), c(n), iQP, jQP, iRPDiv, jRPDiv, weights, order2, singularityMaskChunk, IDMaskChunk);
                Z(iChunk, jChunk) = ZBlock;
            end

            Z = 1i * Z0 / (4 * pi * k) .* Z .* jacobians;
        end
    end

    methods (Access = private)
        function ZBlock = computeZBlock(obj, strategyContext, r, c, iQP, jQP, iRPDiv, jRPDiv, W4D, o2, SM, IDM)
            k = strategyContext.constants.k;

            ZBlock4 = zeros(r, c, 4);
            RPDiv = zeros(r, c, 4);
            W = repmat(W4D, r, c);

            for n = 1:strategyContext.constants.polaritiesCount
                p1 = strategyContext.constants.polarities(n, 1);
                p2 = strategyContext.constants.polarities(n, 2);

                iQPp = iQP(:, :, p1);
                jQPp = jQP(:, :, p2);
                iRPDivp = iRPDiv(:, :, p1);
                jRPDivp = jRPDiv(:, :, p2);

                SMp = repelem(SM(:, :, n), o2, o2);

                D = mom.DistMatrix.getChunk(iQPp, jQPp);
                D0 = (D < eps);

                RPDiv(:, :, n) = iRPDivp.' * jRPDivp;

                G = exp(-1i * k * D) ./ D;
                G(SMp) = G(SMp) - 1 ./ D(SMp);

                I = G;
                I(D0) = -1i * k;

                IW = I .* W;
                Ir = reshape(IW, o2, r, o2, c);
                SS = sum(sum(Ir, 1), 3);
                ZBlock4(:, :, n) = reshape(SS, r, c);
            end

            ZBlock4 = (ZBlock4 + obj.computeElSingularities(strategyContext, SM, IDM)) .* RPDiv;
            ZBlock = sum(ZBlock4, 3);
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

        function singularities = evalElSingularities(~, xs)
            singularities = mom.ImpedanceMatrix.singularity.type3(xs(1, :), xs(2, :));
        end
    end
end