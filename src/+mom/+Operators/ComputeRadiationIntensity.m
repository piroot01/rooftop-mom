classdef ComputeRadiationIntensity < handle
    methods(Static)
        function U = computeInDirection(strategyContext, I, k, rHat)
            F = mom.Operators.ComputeFarField.computeInDirection(strategyContext, k, rHat);
            v = F * I;
            U = 0.5 * sum(abs(v).^2);
        end
    end
end
