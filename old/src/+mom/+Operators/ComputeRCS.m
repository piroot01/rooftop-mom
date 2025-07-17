classdef ComputeRCS < handle
    methods(Static)
        function rho = computeFromExcitation(strategyContext, I, V, k, Ei)
            rho = (strategyContext.constants.Z0 * k)^2 / (4 * pi) * abs(V(:)' * I(:))^2 / norm(Ei)^2;
        end
        
        % probably does not work
        function rho = computeInDirection(strategyContext, I, k, rHat, Ei)
            F = mom.Operators.ComputeFarField.computeInDirection(strategyContext, k, rHat);
            v = F * I;
            rho = 4 * pi * sum(abs(v).^2) / norm(Ei)^2;
        end

        %{
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
        %}
    end
end

