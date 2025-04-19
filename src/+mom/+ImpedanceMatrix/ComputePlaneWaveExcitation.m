classdef ComputePlaneWaveExcitation < handle
    methods (Static)
        function V = compute(strategyContext, direction, polarization)
            k = strategyContext.constants.k;

            order = strategyContext.quadrature.getOrder();
            blockSize = order^2;
            basisFunctionCount = strategyContext.basisFunctions.getCount();

            quadratureWeights = strategyContext.quadrature.getWeights2D();
            quadraturePoints = strategyContext.rooftopBasisFunctions.getQuadraturePoints();
            rooftopPoints = strategyContext.rooftopBasisFunctions.getRooftopPoints();

            exponent = (1i) * k * pagemtimes(direction, quadraturePoints);

            integralKernel = exp(exponent);

            integrand = rooftopPoints .* integralKernel;

            weights = repmat(quadratureWeights, basisFunctionCount, 1);

            integral = pagemtimes(polarization, (integrand .* weights.'));

            integralReshaped = reshape(integral, blockSize, basisFunctionCount, 2);
            V = sum(sum(integralReshaped, 1), 3);
        end
    end
end
