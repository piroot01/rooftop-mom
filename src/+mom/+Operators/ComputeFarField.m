classdef ComputeFarField < handle
    methods(Static)
        function F = computeInDirection(strategyContext, k, rHat)
            % integral constant
            integralConstant = (-1i * k * sqrt(strategyContext.constants.Z0)) / (4 * pi);

            % prepare the phiHat and thetaHat
            zHat = [0; 0; 1];

            rHat = rHat(:);

            phiVec = cross(zHat, rHat);
            phiVecNorm = norm(phiVec);

            if phiVecNorm < eps
                phiHat = [0; 1; 0];
            else
                phiHat = phiVec / phiVecNorm;
            end


            thetaHat = cross(phiHat, rHat);

            % compute the integral
            rooftops = strategyContext.rooftopBasisFunctions;
            quadrature = strategyContext.quadrature;

            order = quadrature.getOrder();
            N = strategyContext.basisFunctions.getCount();
            w2D = quadrature.getWeights2D();
            weights = repmat(w2D, N, 1);

            allQuadraturePoints = rooftops.getQuadraturePoints(false);
            allRooftopPoints = rooftops.getRooftopPoints();

            exponent = pagemtimes(rHat.', allQuadraturePoints);
            phase = exp(1i * k * exponent);

            integrand = allRooftopPoints .* phase;

            integralTheta = pagemtimes(thetaHat.', (integrand .* weights.'));
            integralPhi = pagemtimes(phiHat.', (integrand .* weights.'));

            integralThetaReshaped = reshape(integralTheta, order^2, N, 2);
            integralPhiReshaped = reshape(integralPhi, order^2, N, 2);

            FTheta = sum(sum(integralThetaReshaped, 1), 3);
            FPhi = sum(sum(integralPhiReshaped, 1), 3);

            FTheta = integralConstant * FTheta;
            FPhi = integralConstant * FPhi;

            % combine into a matrix 2 x N
            F = cat(1, FTheta, FPhi);
        end
    end
end
