function result = generateScalarGreensFunction(wavenumber, distanceMatrix)
% generateScalarGreensFunction: Computes the scalar free-space Green's function.
%
%   Note: This function does not include the common 1/(4*pi) normalization factor.
%
%   Inputs:
%       wavenumber: The scalar wavenumber $k$ of the system.
%       distanceMatrix: An MxN matrix of Euclidean distances $R_{ij}$.
%
%   Outputs:
%       result: An MxN matrix representing the scalar Green's function.

    % Evaluate the scalar Green's function G(R) = exp(-ikR) / R.
    result = exp(-1i * wavenumber * distanceMatrix) ./ distanceMatrix;
end
