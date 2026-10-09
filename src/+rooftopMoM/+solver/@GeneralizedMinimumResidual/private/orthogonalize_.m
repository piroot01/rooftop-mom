function [vector, coefficients] = orthogonalize_(vector, basis, columnCount)
% orthogonalize_: Returns the vector orthogonalized against the first columns of the basis by classical Gram-Schmidt
%   run twice; the second pass removes what rounding left after the first. The products take the whole basis, as a
%   block of its columns would be a copy, and the coefficients of the other columns are set to zero.
%
%   Inputs:
%       vector: Nx1 vector.
%       basis: NxK matrix whose first columnCount columns are orthonormal; the others are ignored.
%       columnCount: Number of the columns to orthogonalize against.
%
%   Outputs:
%       vector: Nx1 vector orthogonal to the first columnCount columns.
%       coefficients: Kx1 coefficients of the removed components, zero beyond columnCount.

    coefficients = basis' * vector;
    coefficients((columnCount + 1):end) = 0;
    vector = vector - basis * coefficients;

    correction = basis' * vector;
    correction((columnCount + 1):end) = 0;
    vector = vector - basis * correction;

    coefficients = coefficients + correction;
end
