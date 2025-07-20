function [Y, activeIndices] = rank1AdditionInverseMatrix(Y, Z, activeIndices, n)
% rank1AdditionInverseMatrix: Performs a rank-1 update of the inverse matrix Y when adding the n-th degree of freedom (DOF).
%
%   Inputs:
%       Y: Perturbed inverse matrix (Z^(-1)) with zeros in inactive rows/cols.
%       Z: Original full matrix of all DOFs.
%       activeIndices: 1×K sorted array of active DOF indices in Y.
%       n: Index of the DOF to add to Y.
%
%   Outputs:
%       Y: Rank-1 updated inverse matrix including the n-th DOF.
%       activeIndices: Updated sorted array of active DOF indices (size K+1).

    zn = Z(:, n);
    xn = Y * zn;
    cn = 1 / (Z(n, n) - zn.' * xn);

    % Update existing active submatrix.
    Y = Y + cn * (xn * xn.');

    % Fill in new row and column corresponding to the n-th DOF.
    Y(:, n) = -cn * xn;
    Y(n, :) = Y(:, n).';
    Y(n, n) = cn;

    % Include the n-th DOF in the activeIndices list.
    activeIndices = insertSorted(activeIndices, n);
end
