function [distanceMatrix, functionMatrix] = generateScalarGreensFunction(wavenumber, pointsA, pointsB)
% generateScalarGreensFunction: Computes the scalar free-space Green's function.
%
%   Note: This function does not include the common 1/(4*pi) normalization factor.
%
%   Inputs:
%       wavenumber: The scalar wavenumber $k$ of the system.
%       pointsA:    A 3xM matrix where each column is a 3D source point.
%       pointsB:    A 3xN matrix where each column is a 3D observation point.
%
%   Outputs:
%       distanceMatrix: An MxN matrix of Euclidean distances $R_{ij}$.
%       functionMatrix: An MxN matrix representing the scalar Green's function.

    distanceMatrix = rooftopMoM.core.impedanceMatrix.utils.generateDistanceMatrix(pointsA, pointsB);

    % Evaluate the scalar Green's function G(R) = exp(-ikR) / R.
    functionMatrix = exp(-1i * wavenumber * distanceMatrix) ./ distanceMatrix;
end
