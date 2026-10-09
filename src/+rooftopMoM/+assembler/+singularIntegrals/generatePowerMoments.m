function moments = generatePowerMoments(sizeX, sizeY, power)
% generatePowerMoments: Returns the seven moments of the kernel R^power over the cell pairs at the touching offsets
%   (0, 0), (1, 0), (0, 1) and (1, 1), in closed form (the near-field note in docs/nearfield, Touching moments). The
%   four-dimensional pair integral reduces to an integral over [-1, 1]^2 of the correlation weights of the rooftop
%   profiles times R^power. On each quadrant of the square the weights are polynomials and R^power is singular at most
%   at a corner, so the integral is a sum of integrals of u^m v^n R^power over rectangles of the first quadrant.
%
%   Inputs:
%       sizeX, sizeY: Cell sizes a and b.
%       power: -1, 1 or 3.
%
%   Output:
%       moments: 7x4 moments divided by (a * b)^2, the columns the offsets, the rows M00, M10, M01, M11, N10, N01 and
%           N11: M_alpha_beta weights the pair by xi^alpha xi'^beta and N_alpha_beta by eta^alpha eta'^beta, with the
%           local coordinates xi, eta of the cell at the origin and xi', eta' of the other cell (Cell-pair moments).

    % The y moments are the x moments of the transposed cells at the transposed offsets.
    alongX = generateAlongX_(sizeX, sizeY, power);
    alongY = generateAlongX_(sizeY, sizeX, power);
    moments = [alongX; alongY(2:4, [1, 3, 2, 4])];
end

function moments = generateAlongX_(a, b, power)
% generateAlongX_: Returns the moments M00, M10, M01 and M11 at the four offsets as a 4x4 matrix. They are a^power
%   times those of the cell 1 x ratio, ratio = b / a (Homogeneity).

    ratio = b / a;
    rectangles = integrateRectangles_(power, ratio);

    % Correlation weights W00, W10, W01 and W11 of the profiles (Weights of the rooftop moments), ascending coefficients
    % in nu on the parts [-1, 0] and [0, 1].
    weights = {
        [1, 1], [1, -1]
        [1/2, 1, 1/2], [1/2, 0, -1/2]
        [1/2, 0, -1/2], [1/2, -1, 1/2]
        [1/3, 1/2, 0, -1/6], [1/3, -1/2, 0, 1/6]
    };
    offsets = [0, 1, 0, 1; 0, 0, 1, 1];
    moments = zeros(4, 4);

    % Every moment sums over the four pieces s x t, with the weight k along u and W00 along v.
    for k = 1:4
        for d = 1:4
            for s = 1:2
                [coefficientsU, intervalU] = composePart_(weights{k, s}, s, offsets(1, d), 1, 4);

                for t = 1:2
                    [coefficientsV, intervalV] = composePart_(weights{1, t}, t, offsets(2, d), ratio, 2);
                    moments(k, d) = moments(k, d) ...
                        + coefficientsU * rectangles(:, :, intervalU, intervalV) * coefficientsV.';
                end
            end
        end
    end

    moments = a^power * moments / ratio;
end

function rectangles = integrateRectangles_(power, ratio)
% integrateRectangles_: Returns the integrals of u^m v^n R^power over [0, 1] or [1, 2] times [0, ratio] or
%   [ratio, 2 ratio], indexed by m, n, the u interval and the v interval. They are the corner differences of
%   integrateMonomial, except on the rectangles at least their longer side away from the origin: there the corner
%   differences of an elongated cell cancel, while the integrand is smooth and Gauss-Legendre integrates it to rounding.

    [U, V] = ndgrid([0, 1, 2], ratio * [0, 1, 2]);
    rectangles = zeros(4, 2, 2, 2);

    for m = 0:3
        for n = 0:1
            corners = rooftopMoM.assembler.singularIntegrals.utils.integrateMonomial(power, m, n, U, V);
            rectangles(m + 1, n + 1, :, :) = diff(diff(corners, 1, 1), 1, 2);
        end
    end

    % The singular points of R^power lie at least the longer side away, so 20 points per direction suffice.
    [points, pointWeights] = rooftopMoM.basisFunctions.utils.generateGaussLegendre1D(20, 0, 1);

    for i = 1:2
        for j = 1:2
            if hypot(i - 1, ratio * (j - 1)) < max(1, ratio)
                continue;
            end

            u = i - 1 + points;
            v = ratio * (j - 1 + points.');
            kernel = sqrt(u.^2 + v.^2).^power;

            for m = 0:3
                for n = 0:1
                    rectangles(m + 1, n + 1, i, j) = ratio * pointWeights.' * (u.^m .* v.^n .* kernel) * pointWeights;
                end
            end
        end
    end
end

function [coefficients, interval] = composePart_(weight, part, offset, cellSize, count)
% composePart_: Returns the weight on one part of [-1, 1] as a polynomial in the distance to the singular point,
%   u = sigma * cellSize * (nu - offset) >= 0, and the interval of u: 1 for [0, cellSize], 2 for
%   [cellSize, 2 cellSize] (The polynomial pieces). The coefficients are ascending and padded to count.

    ends = [-1, 0; 0, 1];
    shifted = ends(part, :) - offset;
    sigma = 1 - 2 * (shifted(2) <= 0);
    interval = min(sigma * shifted) + 1;

    % weight(nu) with nu = sigma * u / cellSize + offset; term holds the coefficients of nu^(j - 1) in u.
    coefficients = zeros(1, count);
    term = 1;

    for j = 1:numel(weight)
        coefficients(1:j) = coefficients(1:j) + weight(j) * term;
        term = conv(term, [offset, sigma / cellSize]);
    end
end
