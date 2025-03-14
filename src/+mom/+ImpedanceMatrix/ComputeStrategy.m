classdef (Abstract) ComputeStrategy < handle
    methods (Abstract)
        Z = compute(obj, strategyContext);
    end

    methods
        function mask = getSingularityMask(~, strategyContext)
            basisFuncsCount = strategyContext.basisFunctions.getCount();
            basisFuncs = strategyContext.basisFunctions.getFaces();

            mask = false(basisFuncsCount);

            for n = 1:(strategyContext.constants.polaritiesCount)
                p1 = strategyContext.constants.polarities(n, 1);
                p2 = strategyContext.constants.polarities(n, 2);

                mask(:, :, n) = (basisFuncs(p1, :).' == basisFuncs(p2, :));
            end
        end

        function mask = getIDMask(~, strategyContext)
            basisFuncIDs = double(strategyContext.basisFunctions.getFuncIDs());
            mask = (basisFuncIDs == basisFuncIDs.') .* (basisFuncIDs.' * ones(1, numel(basisFuncIDs)));
        end

        function jacobians = getJacobians(~, strategyContext)
            basisFuncsIDs = strategyContext.basisFunctions.getFuncIDs();
            IDCount = strategyContext.basisFunctions.getIDCount();

            uniqueJacobians = zeros(1, IDCount);
            for n = 1:IDCount
                [a, b] = strategyContext.basisFunctions.getFaceSize(n);
                uniqueJacobians(n) = a * b;
            end

            jacobians2D = uniqueJacobians(basisFuncsIDs);
            jacobians = jacobians2D.' * jacobians2D;
        end
    end
end