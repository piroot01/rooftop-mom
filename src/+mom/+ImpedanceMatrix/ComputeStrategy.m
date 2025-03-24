classdef (Abstract) ComputeStrategy < handle
    methods (Abstract)
        Z = compute(obj, strategyContext);
    end

    methods
        function mask = getSingularityMask(~, strategyContext)
            p1 = strategyContext.constants.polarities(:, 1);
            p2 = strategyContext.constants.polarities(:, 2);
            polaritiesCount = strategyContext.constants.polaritiesCount;

            basisFunctionCount = strategyContext.basisFunctions.getCount();
            basisFunctions = strategyContext.basisFunctions.getFaces()

            mask = (reshape(basisFunctions(p1, :).', basisFunctionCount, 1, polaritiesCount) == reshape(basisFunctions(p2, :).', 1, basisFunctionCount, polaritiesCount));
        end

        function mask = getIDMask(~, strategyContext)
            basisFuncIDs = double(strategyContext.basisFunctions.getIDs());
            mask = (basisFuncIDs == basisFuncIDs.') .* (basisFuncIDs.' * ones(1, numel(basisFuncIDs)));
        end

        function jacobians = getJacobians(~, strategyContext)
            basisFuncsIDs = strategyContext.basisFunctions.getIDs();
            rooftopUniqueUniqueSizes = strategyContext.rooftopBasisFunctions.getUniqueUniqueRooftopSizes();
            uniqueJacobians = prod(rooftopUniqueUniqueSizes, 2);
            jacobians2D = uniqueJacobians(basisFuncsIDs);
            jacobians = jacobians2D(:) * jacobians2D(:).';
        end
    end
end
