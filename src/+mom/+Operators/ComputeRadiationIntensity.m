classdef ComputeRadiationIntensity < handle
    methods(Static)
        function U = computeFromExcitation(strategyContext, I, V, k)
            v = V(:)' * I(:);
            U = 0.5 * strategyContext.constants.Z0 * k^2 / (4 * pi)^2 * abs(v)^2;
        end

        function U = computeInDirection(strategyContext, I, k, rHat)
            F = mom.Operators.ComputeFarField.computeInDirection(strategyContext, k, rHat);
            v = F * I;
            U = 0.5 * sum(abs(v).^2);
        end
    end
end
