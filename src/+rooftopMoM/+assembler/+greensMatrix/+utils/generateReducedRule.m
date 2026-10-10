function rule = generateReducedRule(order, partCounts)
% generateReducedRule: Returns the quadrature rule of the seven cell-pair moments. The kernel depends on the two cells
%   only through the differences nu = xi - xi' and eta - eta' of their local coordinates, so every moment is an
%   integral over [-1, 1]^2 of the kernel times correlation weights of the rooftop profiles, which are polynomials on
%   [-1, 0] and on [0, 1] (the near-field note in docs/nearfield, Two-dimensional representation). The rule is
%   Gauss-Legendre of the given order on each of partCounts equal parts of [-1, 0] and of [0, 1]:
%   (2 * order)^2 * prod(partCounts) kernel values per cell pair instead of the order^4 of the tensor rule on both
%   cells.
%
%   Inputs:
%       order: Number of Gauss-Legendre points per part.
%       partCounts: 1x2 numbers of parts of [-1, 0] and of [0, 1] in x and in y.
%
%   Output:
%       rule: Struct with the fields:
%           nodesX: Px1 values of nu, Px = 2 * order * partCounts(1).
%           nodesY: 1xPy values of eta, Py = 2 * order * partCounts(2).
%           weights: 7x(Px * Py) weights of the moments M00, M10, M01, M11, N10, N01 and N11
%               (singularIntegrals.generatePowerMoments), nu running fastest.

    [nodesX, weightsX] = generateAxis_(order, partCounts(1));
    [nodesY, weightsY] = generateAxis_(order, partCounts(2));

    rule.nodesX = nodesX;
    rule.nodesY = nodesY.';

    % The weight of the node (nu_k, eta_l) is the product of a column of weightsX at k and one of the columns of
    % weightsY at l, see Theorem 4.1. M weights nu by W00, W10, W01 or W11 and eta just by W00, N does it the other
    % way round: with the x and y axes exchanged, M becomes N (Lemma 2.1 (iii)). Thus
    %   kron(weightsY(:, 1), weightsX) is W00(eta) times the four columns of nu: M00, M10, M01 and M11,
    %   kron(weightsY(:, 2:4), weightsX(:, 1)) is W00(nu) times the last three columns of eta: N10, N01 and N11.
    % Because M00 and N00 are the same, only one of them is kept. The transpose leaves a single moment in each row.
    rule.weights = [kron(weightsY(:, 1), weightsX), kron(weightsY(:, 2:4), weightsX(:, 1))].';
end

function [nodes, weights] = generateAxis_(order, partCount)
% generateAxis_: Returns the nodes on [-1, 1] and the Gauss-Legendre weights times the correlation weights W00, W10,
%   W01 and W11 of the profiles 1 and xi (Weights of the rooftop moments).

    % points is an order x 1 column on [0, 1]. Adding the row 0:(partCount - 1) copies it onto the intervals [0, 1],
    % ..., [partCount - 1, partCount], an order x partCount matrix. Dividing by partCount maps them back onto [0, 1],
    % now split into partCount equal parts of the right half I+ (Corollary 4.1), and reshape stacks them into one
    % column.
    [points, pointWeights] = rooftopMoM.basisFunctions.utils.generateGaussLegendre1D(order, 0, 1);
    right = reshape((points + (0:(partCount - 1))) / partCount, [], 1);

    % The left half I- = [-1, 0] gets the same parts, shifted by -1. Splitting [-1, 1] at 0 keeps every part inside one
    % half, where the weights are polynomials, so their kink at nu = 0 never falls inside a part. Each part has the
    % length 1 / partCount, which scales its weights, and there are 2 * partCount parts in total.
    nodes = [right - 1; right];
    pointWeights = repmat(pointWeights / partCount, 2 * partCount, 1);

    % W10 and W01 of (23) are not even, W01(nu) = W10(-nu), so they take one polynomial on each half.
    isLeft = nodes < 0;
    weight10 = (1 - nodes.^2) / 2;
    weight10(isLeft) = (1 + nodes(isLeft)).^2 / 2;
    weight01 = (1 - nodes).^2 / 2;
    weight01(isLeft) = (1 - nodes(isLeft).^2) / 2;

    % W00 and W11 are even, so one formula in |nu| covers both halves: W00 = 1 - |nu|, and
    % W11 = (1 - |nu|)^2 (2 + |nu|) / 6, which is (2 - 3 nu + nu^3) / 6 on I+ and (1 + nu)^2 (2 - nu) / 6 on I-.
    distance = abs(nodes);

    % Every column is the Gauss-Legendre weight times one correlation weight, W00, W10, W01 and W11.
    weights = pointWeights .* [1 - distance, weight10, weight01, (1 - distance).^2 .* (2 + distance) / 6];
end
