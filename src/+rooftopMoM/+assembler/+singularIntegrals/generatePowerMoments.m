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
    % in nu on the parts [-1, 0] and [0, 1]. See Corollary 4.1 (Weights of the rooftop moments).
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
                % The weight k on the part s along u: the coefficients c^s_m of (54), i.e. one cell of Table 2, and
                % the interval A_s, 1 for [0, a] or 2 for [a, 2a].
                [coefficientsA, intervalA] = composePart_(weights{k, s}, s, offsets(1, d), 1, 4);

                for t = 1:2
                    % W00 on the part t along v: the coefficients e^t_n of (54) and the interval B_t, 1 for [0, b] or
                    % 2 for [b, 2b].
                    [coefficientsB, intervalB] = composePart_(weights{1, t}, t, offsets(2, d), ratio, 2);

                    % The piece s x t of (55), the sum over m and n of c^s_m e^t_n times the integral of
                    % u^m v^n R^power over A_s x B_t: the row c^s (1x4) times these integrals (4x2) times the column
                    % e^t (2x1).
                    moments(k, d) = moments(k, d) ...
                        + coefficientsA * rectangles(:, :, intervalA, intervalB) * coefficientsB.';
                end
            end
        end
    end

    % The factor 1 / (ab) of (55) is 1 / ratio for the cell 1 x ratio, and a^power restores the size (Corollary 6.1).
    moments = a^power * moments / ratio;
end

function rectangles = integrateRectangles_(power, ratio)
% integrateRectangles_: Returns the integrals of u^m v^n R^power over [0, 1] or [1, 2] times [0, ratio] or
%   [ratio, 2 ratio], indexed by m, n, the u interval and the v interval. They are the corner differences of
%   integrateMonomial, except on the rectangles at least their longer side away from the origin: there the corner
%   differences of an elongated cell cancel, while the integrand is smooth and Gauss-Legendre integrates it to rounding.

    % Start by gathering the 9 corners of the four rectangles A_s x B_t, see Figure 6. By the homogeneity of M, see
    % Corollary 6.1, a is factored out, so the cell is 1 x ratio with ratio = b / a.
    [as, bs] = ndgrid([0, 1, 2], ratio * [0, 1, 2]);

    % The first dimension holds the 4 powers of u (0, 1, 2, 3), the second the 2 powers of v (0, 1), the third the
    % interval A_s along u (1 <=> [0, a], 2 <=> [a, 2a]) and the last one the interval B_t along v (1 <=> [0, b],
    % 2 <=> [b, 2b]).
    rectangles = zeros(4, 2, 2, 2);

    for m = 0:3
        for n = 0:1
            corners = rooftopMoM.assembler.singularIntegrals.utils.integrateMonomial(power, m, n, as, bs);

            % corners is a 3x3 matrix holding Gmn(as(k, l), bs(k, l)) at the position (k, l). To obtain the four
            % results of (56), the first diff takes Gmn(as(k + 1, l), ...) - Gmn(as(k, l), ...), giving a 2x3 matrix,
            % and the second diff takes Gmn(..., bs(k, l + 1)) - Gmn(..., bs(k, l)), giving the final 2x2 matrix.
            rectangles(m + 1, n + 1, :, :) = diff(diff(corners, 1, 1), 1, 2);
        end
    end

    % The closed form is exact, but on the rectangles away from the origin it is a difference of large, nearly
    % cancelling terms (Remark after Theorem 6.1) and loses about the square of the aspect ratio in relative accuracy,
    % e.g. 1e-9 instead of 1e-15 at b / a = 1e-3. There the integrand is smooth and the singular point lies at least
    % the longer side away, so Gauss-Legendre with 20 points per direction integrates it to rounding.
    [points, pointWeights] = rooftopMoM.basisFunctions.utils.generateGaussLegendre1D(20, 0, 1);

    for i = 1:2
        for j = 1:2
            % Check whether the rectangle is at least its longer side away from the origin. Only then is it far enough
            % from the singularity to apply Gauss-Legendre safely; otherwise it keeps the closed form.
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
%
%   Each call is one cell of Table 2, computed numerically instead of stored: weight and part pick the row and the
%   piece, offset the column, and the result is the table entry expanded in u with a = cellSize (1 along u, ratio
%   along v). sigma and interval are the sigma and u rows of the table header.

    % The ends of the part [-1, 0] or [0, 1], shifted so that the singular value nu = offset lies at 0.
    ends = [-1, 0; 0, 1];
    shifted = ends(part, :) - offset;

    % On the part, nu - offset keeps one sign: sigma = -1 if the shifted part ends at or below 0, +1 otherwise, so
    % that sigma * (nu - offset) = |nu - offset| >= 0, see (52).
    sigma = 1 - 2 * (shifted(2) <= 0);

    % sigma * shifted holds the ends of |nu - offset| = u / cellSize, either [0, 1] or [1, 2], see (53). Its lower
    % end, 0 or 1, gives the interval index 1 or 2.
    interval = min(sigma * shifted) + 1;

    % Substitute nu = offset + sigma * u / cellSize into weight(nu) = sum_j weight(j) nu^(j - 1). term holds the
    % ascending coefficients in u of nu^(j - 1): it starts at nu^0 = 1, and every step multiplies it by
    % nu = [offset, sigma / cellSize] through conv, the product of two polynomials, raising its degree by one. At
    % step j it has j coefficients, and weight(j) times it is added to the first j of the result. E.g. W00 on [-1, 0]
    % at offset 1 is 1 + nu with nu = 1 - u, which gives [2, -1], i.e. 2 - u. The coefficients are padded to the
    % polynomial size with zeros.
    coefficients = zeros(1, count);
    term = 1;

    for j = 1:numel(weight)
        coefficients(1:j) = coefficients(1:j) + weight(j) * term;
        term = conv(term, [offset, sigma / cellSize]);
    end
end
