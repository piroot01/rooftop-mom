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
        function Z = compute(obj, strategyContext, gridSize)
            computeStrategyOptimiser = mom.ImpedanceMatrix.ComputeStrategyOptimiser(strategyContext.basisFunctions, gridSize);
            [mask, interactions] = computeStrategyOptimiser.getMaskAndInteractions();
            chunks = obj.computeChunkedContext_.generateChunksSpecificInteractions(strategyContext.quadrature.getOrder(), interactions, strategyContext.basisFunctions.getCount());

            ZUniqueValues = obj.computeStrategy_.compute(strategyContext, chunks);
            Z = computeStrategyOptimiser.applyMask(ZUniqueValues, mask, interactions);
        end
    end
end
