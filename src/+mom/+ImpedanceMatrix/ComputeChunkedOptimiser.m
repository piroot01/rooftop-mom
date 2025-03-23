classdef ComputeChunkedOptimiser < handle
    properties
        computeStrategy_;
        computeChunkedContext_ mom.ImpedanceMatrix.ComputeChunkedContext;
    end

    methods
        function obj = ComputeChunkedOptimiser(computeStrategy, computeChunkedContext)
            obj.computeStrategy_ = computeStrategy;
            obj.computeChunkedContext_ = computeChunkedContext;
        end
    end

    methods
        % NOTE: does not work with multiple antennas yet
        function Z = compute(obj, strategyContext, gridSizes)
            computeStrategyOptimiser = mom.ImpedanceMatrix.ComputeStrategyOptimiser(strategyContext.basisFunctions, strategyContext.rooftopBasisFunctions.getUniqueUniqueRooftopSizes(), gridSizes);
            [masksCell, interactionsCell] = computeStrategyOptimiser.getAllMasksAndInteractions();

            [mask, interactions, chunks] = obj.generateMaskAndInteractionsAndChunks(strategyContext, masksCell, interactionsCell);

            celldisp(chunks);

            ZUniqueValues = obj.computeStrategy_.compute(strategyContext, chunks);

            Z = computeStrategyOptimiser.applyMask(ZUniqueValues, mask, interactions);
        end
    end

    methods (Access = private)
        function [mask, interactions, chunks] = generateMaskAndInteractionsAndChunks(obj, strategyContext, masksCell, interactionsCell)
            basisFunctions = strategyContext.basisFunctions.getFacesAsCell();
            IDCount = strategyContext.basisFunctions.getIDCount();
            IDMask = obj.computeStrategy_.getIDMask(strategyContext);
            quadratureOrder = strategyContext.quadrature.getOrder();

            mask = zeros(size(IDMask), 'uint32');
            interactions = [interactionsCell{:}];

            for n = 1:IDCount
                basisFunctionCount = size(basisFunctions{n}, 2);
                chunks = obj.computeChunkedContext_.generateChunksSpecificInteractions(quadratureOrder, interactionsCell{n}, basisFunctionCount);
                mask(IDMask == n) = masksCell{n};
            end

            mask
        end
    end
end
