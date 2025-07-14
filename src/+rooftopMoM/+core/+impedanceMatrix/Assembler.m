classdef Assembler < handle
    methods (Abstract)
        Z = generate(obj, assemblerContext);
    end

    methods
        function mask = generateSingularityMask(~, assemblerContext)
            p1 = [1, 1, 2, 2];
            p2 = [1, 2, 1, 2];
            pCount = 4;

            masks = cellfun( ...
                @(idx) (reshape(idx(p1, :).', basisSize, 1, size(idx, 2), pCount) == reshape(idx(p2, :).', 1, size(idx, 2), pCount)), ...
                assemblerContext.indices, 'UniformOutput', false)

            mask = blkdiag(masks)
        end
    end
end


%{
        function mask = getSingularityMask(~, strategyContext)
            p1 = strategyContext.constants.polarities(:, 1);
            p2 = strategyContext.constants.polarities(:, 2);
            polaritiesCount = strategyContext.constants.polaritiesCount;

            basisFunctionCount = strategyContext.basisFunctions.getCount();
            basisFunctions = strategyContext.basisFunctions.getFaces();

            mask = (reshape(basisFunctions(p1, :).', basisFunctionCount, 1, polaritiesCount) == reshape(basisFunctions(p2, :).', 1, basisFunctionCount, polaritiesCount));
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
%}
