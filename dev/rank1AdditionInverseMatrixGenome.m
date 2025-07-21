function [Y, genome] = rank1AdditionInverseMatrixGenome(Y, Z, genome, n)
% rank1AdditionInverseMatrixGenome: Performs a rank-1 update of the inverse matrix Y when adding the n-th degree of freedom (DOF).
%
%   Inputs:
%       Y: Perturbed inverse matrix (Z^(-1)) with zeros in inactive rows/cols.
%       Z: Original full matrix of all DOFs.
%       genome: 1×N logical array of active DOFs.
%       n: Index of the DOF to add to Y.
%
%   Outputs:
%       Y: Rank-1 updated inverse matrix including the n-th DOF.
%       genome: Updated genome of active DOF indices.

    zn = Z(:, n);
    xn = Y * zn;
    cn = 1 / (Z(n, n) - zn.' * xn);

    % Update existing active submatrix.
    Y = Y + cn * (xn * xn.');

    % Fill in new row and column corresponding to the n-th DOF.
    Y(:, n) = -cn * xn;
    Y(n, :) = Y(:, n).';
    Y(n, n) = cn;

    % Activate the n-th DOF.
    genome(n) = true;
end
