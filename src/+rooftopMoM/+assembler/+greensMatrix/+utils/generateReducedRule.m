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

    % M weights nu by Wij and eta by W00, N the other way round.
    rule.nodesX = nodesX;
    rule.nodesY = nodesY.';
    rule.weights = [kron(weightsY(:, 1), weightsX), kron(weightsY(:, 2:4), weightsX(:, 1))].';
end

function [nodes, weights] = generateAxis_(order, partCount)
% generateAxis_: Returns the nodes on [-1, 1] and the Gauss-Legendre weights times the correlation weights W00, W10,
%   W01 and W11 of the profiles 1 and xi (Weights of the rooftop moments).

    [points, pointWeights] = rooftopMoM.basisFunctions.utils.generateGaussLegendre1D(order, 0, 1);
    right = reshape((points + (0:(partCount - 1))) / partCount, [], 1);
    nodes = [right - 1; right];
    pointWeights = repmat(pointWeights / partCount, 2 * partCount, 1);

    isLeft = nodes < 0;
    weight10 = (1 - nodes.^2) / 2;
    weight10(isLeft) = (1 + nodes(isLeft)).^2 / 2;
    weight01 = (1 - nodes).^2 / 2;
    weight01(isLeft) = (1 - nodes(isLeft).^2) / 2;
    distance = abs(nodes);

    weights = pointWeights .* [1 - distance, weight10, weight01, (1 - distance).^2 .* (2 + distance) / 6];
end
