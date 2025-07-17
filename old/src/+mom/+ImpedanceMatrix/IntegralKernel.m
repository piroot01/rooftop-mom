classdef IntegralKernel < handle
    methods (Static)
        function integralKernel = compute(strategyContext, rowQuadraturePointsPolarity, columnQuadraturePointsPolarity)
            k = strategyContext.constants.k;

            integralKernel.distanceMatrix = mom.DistanceMatrix.compute(rowQuadraturePointsPolarity, columnQuadraturePointsPolarity);
            integralKernel.greensFunction = exp((-1i) * k * integralKernel.distanceMatrix) ./ integralKernel.distanceMatrix;
        end
    end
end
