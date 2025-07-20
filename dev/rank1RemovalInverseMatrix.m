function [Y, activeIndices] = rank1RemovalInverseMatrix(Y, activeIndices, n)
% rank1RemovalInverseMatrix: Performs a rank-1 update of the inverse matrix Y when removing the n-th degree of freedom (DOF).
%
%   Inputs:
%       Y: Perturbed inverse matrix (Z^(-1)) with zeros in inactive rows/cols.
%       activeIndices: 1×K sorted array of active DOF indices in Y.
%       n: Index of the DOF to remove from Y.
%
%   Outputs:
%       Y: Rank-1 updated inverse matrix excluding the n-th DOF.
%       activeIndices: Updated sorted array of active DOF indices (size K-1).

    yn = Y(:, n);

    % Perform the basis function removal.
    Y = Y - (yn / Y(n, n)) * yn.';

    % Update the indices.
    activeIndices(n) = [];
end
