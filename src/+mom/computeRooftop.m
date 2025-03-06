function [rooftopP, rooftopN] = computeRooftop(basisFunctionMesh, quadraturePoints)
    % input:
    %   basisFunctionMesh: 3x8 matrix (8 vertices - technically 6 is enough)
    %   polarity: 0 = positive, 1 = negative
    %   quadraturePoints: 3xN matrix (N is the order of the qudrature squared)

    v1p = basisFunctionMesh(:, 1);
    v1n = basisFunctionMesh(:, 5);
    v2n = basisFunctionMesh(:, 6);

    % Compute the difference vector and its norm.
    d1 = v1n - v1p;
    d2 = v2n - v1n;
    l1 = sqrt(sum(d1 .* d1));
    l2 = sqrt(sum(d2 .* d2));
    e1 = d1 / l1;
    e2 = d2 / l2;

    % Compute the rooftop using a vectorized projection.
    rooftopP = (1 / l1) * (e1 * (e1.' * (quadraturePoints(:, 1:4) - v1p)));
    rooftopN = (-1 / l2) * (e2 * (e2.' * (quadraturePoints(:, 5:8) - v2n)));
end