classdef ComputeRCS < handle
    methods(Static)
        function [RCS_dB, thetaVec, phiVec] = compute(strategyContext, I, f, thetaSamples, phiSamples)
            if nargin < 3
                f = 0;
            end

            if nargin < 4 || isempty(thetaSamples)
                thetaSamples = 181;
            end

            if nargin < 5 || isempty(phiSamples)
                phiSamples   = 361;
            end

            c0 = 299792458;

            lambda = c0 / f;

            k      = 2 * pi / lambda;
            rbf    = strategyContext.rooftopBasisFunctions;
            quad   = strategyContext.quadrature;

            order  = quad.getOrder();
            N      = strategyContext.basisFunctions.getCount();
            w2D    = quad.getWeights2D();
            weightsAll = repmat(w2D, N, 1);

            quadPtsAll = rbf.getQuadraturePoints();
            rrAll      = rbf.getRooftopPoints();

            thetaVec = linspace(0, pi, thetaSamples);
            phiVec   = linspace(0, 2*pi, phiSamples);

            RCS_lin = zeros(phiSamples, thetaSamples);

            for p = 1:phiSamples
                phi = phiVec(p);
                for t = 1:thetaSamples
                    th = thetaVec(t);

                    hatr = [sin(th)*cos(phi), sin(th)*sin(phi), cos(th)];  

                    exponent = pagemtimes(hatr, quadPtsAll);
                    phase = exp(-1i * k * exponent);  

                    integrandVec = rrAll .* phase;
                    integralKernel = integrandVec .* weightsAll.';

                    integralKernel1 = integralKernel(1, :, :);
                    integralKernel2 = integralKernel(2, :, :);
                    integralKernel3 = integralKernel(3, :, :);

                    integrals1 = reshape(integralKernel1, order^2, N, 2);
                    integrals2 = reshape(integralKernel2, order^2, N, 2);
                    integrals3 = reshape(integralKernel3, order^2, N, 2);

                    integralsSummed1 = sum(sum(integrals1, 1), 3);
                    integralsSummed2 = sum(sum(integrals2, 1), 3);
                    integralsSummed3 = sum(sum(integrals3, 1), 3);

                    integralsSummed = [integralsSummed1; integralsSummed2; integralsSummed3];

                    E_sca = sum(I(:).' .* integralsSummed, 2);

                    %— Compute monostatic RCS (E_inc=1):
                    RCS_lin(p,t) = 4*pi * norm(E_sca).^2;
                end
            end

            RCS_dB = 10*log10(RCS_lin + eps);
        end

        function U = computeU(strategyContext, I, k, theta, phi)
            rbf    = strategyContext.rooftopBasisFunctions;
            quad   = strategyContext.quadrature;

            order  = quad.getOrder();
            N      = strategyContext.basisFunctions.getCount();
            w2D    = quad.getWeights2D();
            weightsAll = repmat(w2D, N, 1);

            quadPtsAll = rbf.getQuadraturePoints(false);
            rrAll      = rbf.getRooftopPoints();

            hatr = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]; 

            exponent = pagemtimes(hatr, quadPtsAll);
            phase = exp(1i * k * exponent);

            integrand = rrAll .* phase;

            thetahat = [cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta)];
            phihat   = [-sin(phi), cos(phi), 0];

            integralTheta = pagemtimes(thetahat, (integrand .* weightsAll.'));
            integralPhi = pagemtimes(phihat, (integrand .* weightsAll.'));

            integralThetaReshaped = reshape(integralTheta, order^2, N, 2);
            FTheta = sum(sum(integralThetaReshaped, 1), 3);
            FTheta = (-1i * strategyContext.constants.Z0 * k / (4 * pi)) * FTheta;

            integralPhiReshaped = reshape(integralPhi, order^2, N, 2);
            FPhi = sum(sum(integralPhiReshaped, 1), 3);
            FPhi = (-1i * strategyContext.constants.Z0 * k / (4 * pi)) * FPhi;

            UTheta = (1/(2 * strategyContext.constants.Z0)) * abs(FTheta * I)^2;
            UPhi = (1/(2 * strategyContext.constants.Z0)) * abs(FPhi * I)^2;

            U = UTheta + UPhi;
        end
    end
end

