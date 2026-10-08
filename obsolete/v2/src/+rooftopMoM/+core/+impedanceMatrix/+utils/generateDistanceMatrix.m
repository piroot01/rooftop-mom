function result = generateDistanceMatrix(pointsA, pointsB)
% generateDistanceMatrix: Computes the Euclidean distance matrix between two sets of points.
%
%   Inputs:
%       pointsA: A 3xMx2 matrix where each column represents a 3D point.
%       pointsB: A 3xNx2 matrix where each column represents a 3D point.
%
%   Output:
%       result: An MxN matrix where `result(i, j)` is the Euclidean distance between `pointsA(:, i, :)` and `pointsB(:, j, :)`.

    % Calculate the squared L2 norm for each point.
    normsSqA = sum(pointsA.^2, 1);
    normsSqB = sum(pointsB.^2, 1);

    % Compute all dot products via matrix multiplication.
    dotProduct = pagemtimes(permute(pointsA, [2, 1, 3]), pointsB);

    % Finally use the identity:
    %   $ \| \mathbf{a} - \mathbf{b} \|^2 = \| \mathbf{a} \|^2 + \| \mathbf{b} \|^2 - 2 \mathbf{a} \cdot \mathbf{b} $
    result = sqrt(permute(normsSqA, [2, 1, 3]) + normsSqB - 2 * dotProduct);
end
