function [points, weights] = generateGaussLegendre1D(order, a, b)
% generateGaussLegendre1D: Returns the Gauss-Legendre points and weights on the interval [a, b].
%   The points are the eigenvalues of the Jacobi matrix of the Legendre polynomials, the weights follow from the
%   first components of its eigenvectors (Golub-Welsch algorithm).
%
%   Inputs:
%       order: Positive integer specifying the number of points. Polynomials up to the degree 2 * order - 1 are
%           integrated exactly.
%       a: Start of the interval.
%       b: End of the interval.
%
%   Outputs:
%       points: order x 1 vector of the points in ascending order.
%       weights: order x 1 vector of the weights.

    n = (1:(order - 1)).';
    offDiagonal = n ./ sqrt(4 * n.^2 - 1);
    [vectors, values] = eig(diag(offDiagonal, 1) + diag(offDiagonal, -1), 'vector');
    [values, permutation] = sort(values);

    % The weights on [-1, 1] sum to 2.
    weights = 2 * vectors(1, permutation).'.^2;

    points = (a + b) / 2 + (b - a) / 2 * values;
    weights = (b - a) / 2 * weights;
end
