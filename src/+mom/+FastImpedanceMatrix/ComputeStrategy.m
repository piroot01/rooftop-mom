classdef (Abstract) ComputeStrategy < handle
    methods
        function mask = getSingularityMask(~, strategyContext, direction1, direction2, ID)
            p1 = strategyContext.constants.polarities(:, 1);
            p2 = strategyContext.constants.polarities(:, 2);
            polaritiesCount = strategyContext.constants.polaritiesCount;

            basis1 = strategyContext.basisFunctions.getFacesInDirection(direction1);
            basisFunctionFaces1 = [basis1{ID}];

            basis2 = strategyContext.basisFunctions.getFacesInDirection(direction2);
            basisFunctionFaces2 = [basis2{ID}];

            mask = (reshape(basisFunctionFaces1(p1, :).', size(basisFunctionFaces1, 2), 1, polaritiesCount) == reshape(basisFunctionFaces2(p2, :).', 1, size(basisFunctionFaces2, 2), polaritiesCount));
        end

        function mask = getIDMask(~, strategyContext)
            basisFuncIDs = double(strategyContext.basisFunctions.getIDs());
            mask = (basisFuncIDs == basisFuncIDs.') .* (basisFuncIDs.' * ones(1, numel(basisFuncIDs)));
        end

        function jacobians = getJacobians(~, strategyContext)
            rooftopSizes = strategyContext.rooftopBasisFunctions.getRooftopSizes();
            jacobians2D = prod(rooftopSizes, 2);
            jacobians = jacobians2D(:) * jacobians2D(:).';
        end
    end
end
