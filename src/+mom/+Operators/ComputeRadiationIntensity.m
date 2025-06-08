classdef ComputeRadiationIntensity < handle
    methods(Static)
        function U = computeInDirection(strategyContext, I, k, rHat)
            F = mom.Operators.ComputeFarField.computeInDirection(strategyContext, k, rHat);
            U = 1/2 * real((I' * F') * (F * I));
        end
    end
end
