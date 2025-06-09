classdef ComputeSphericalProjector < handle
    methods(Static)
        function S = compute(strategyContext, k, L)
            % integral constant
            integralConstant = k * sqrt(strategyContext.constants.Z0);

            % compute the integral
            rooftops = strategyContext.rooftopBasisFunctions;
            quadrature = strategyContext.quadrature;

            order = quadrature.getOrder();
            N = strategyContext.basisFunctions.getCount();
            w2D = quadrature.getWeights2D();
            weights = repmat(w2D, N, 1);

            allQuadraturePoints = rooftops.getQuadraturePoints(false);
            allRooftopPoints = rooftops.getRooftopPoints();

            % convert to spherical coordinates
            x = allQuadraturePoints(1, :, :);
            y = allQuadraturePoints(2, :, :);
            z = allQuadraturePoints(2, :, :);

            r = sqrt(x.^2 + y.^2 + z.^2);
            theta = acos(z ./ r);
            phi = atan2(y, x);

            directions = [phi(:), theta(:)];
            allY = mom.Operators.ComputeSphericalProjector.computeRealScalarSphericalHarmonics(L, directions)

            Nalpha = 2 * L * (L + 2);
            S = zeros(Nalpha, N);
            idx = 0;

            for l = 1:L
                for m = 0:l
                    for sigma = ['e', 'o']
                        for tau = 1:2
                            idx = idx + 1;
                        end
                    end
                end
            end

            idx
        end
    end

    methods (Static, Access = private)
        function Y = computeRealScalarSphericalHarmonics(N, directions)
            % credit to polarch's function getSH

            Ndirs = size(directions, 1);
            Nharm = (N + 1)^2;

            Y = zeros(Nharm, Ndirs);
            idx = 0;

            for n = 0:N
                m = (0:n).';

                % vector of unnormalised associated Legendre functions of current order
                Lnm_real = legendre(n, cos(directions(:, 2).'));
                if n ~= 0
                    % cancel the Condon-Shortley phase from the definition of
                    % the Legendre functions to result in signless real SH
                    condon = (-1).^[m(end:(-1):2);m] * ones(1,Ndirs);
                    Lnm_real = condon .* [Lnm_real(end:(-1):2, :); Lnm_real];
                end

                % normalisations
                norm_real = sqrt( (2 * n + 1)*factorial(n - m) ./ (4 * pi * factorial(n + m)));

                % convert to matrix, for direct matrix multiplication with the rest
                Nnm_real = norm_real * ones(1, Ndirs);
                if n~=0
                    Nnm_real = [Nnm_real(end:-1:2, :); Nnm_real];
                end

                CosSin = zeros(2*n+1,Ndirs);
                % zero degree
                CosSin(n + 1,:) = ones(1, size(directions, 1));
                % positive and negative degrees
                if n ~= 0
                    CosSin(m(2:end) + n + 1, :) = sqrt(2) * cos(m(2:end) * directions(:, 1).');
                    CosSin(-m(end:(-1):2) + n + 1, :) = sqrt(2) * sin(m(end:-1:2) * directions(:, 1).');
                end
                Ynm = Nnm_real .* Lnm_real .* CosSin;

                Y((idx + 1):(idx + (2 * n + 1)), :) = Ynm;
                idx = idx + 2 * n + 1;
            end

            Y = Y.';
        end
    end
end
