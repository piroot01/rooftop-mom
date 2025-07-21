function [Y, genome] = rank1RemovalInverseMatrixGenome(Y, genome, n)
% rank1RemovalInverseMatrixGenome: Performs a rank-1 update of the inverse matrix Y when removing the n-th degree of freedom (DOF).
%
%   Inputs:
%       Y: Perturbed inverse matrix (Z^(-1)) with zeros in inactive rows/cols.
%       genome: 1×N logical array of active DOFs.
%       n: Index of the DOF to remove from Y.
%
%   Outputs:
%       Y: Rank-1 updated inverse matrix excluding the n-th DOF.
%       genome: Updated genome of active DOF indices.

    yn = Y(:, n);

    % Perform the basis function removal.
    Y = Y - (yn / Y(n, n)) * yn.';

    % Update the indices.
    genome(n) = false;
end
