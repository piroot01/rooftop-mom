function [forwardX, backwardX, forwardY, backwardY] = generateSpectralShiftOperators(circSizeX, circSizeY)
% generateSpectralShiftOperators: Computes the phase multipliers of unit spatial shifts in the Fourier domain
%   as full matrices. The forward shift is exp(2i * pi * u / circSizeX) (similarly for y), the backward shift
%   is its complex conjugate.
%
%   Inputs:
%       circSizeX: Positive integer specifying the circulant size in the x direction.
%       circSizeY: Positive integer specifying the circulant size in the y direction.
%
%   Outputs:
%       forwardX: circSizeX x circSizeY matrix of the forward shift in the x direction.
%       backwardX: circSizeX x circSizeY matrix of the backward shift in the x direction.
%       forwardY: circSizeX x circSizeY matrix of the forward shift in the y direction.
%       backwardY: circSizeX x circSizeY matrix of the backward shift in the y direction.

    [forwardX, backwardX, forwardY, backwardY] = ...
        rooftopMoM.core.greensMatrix.utils.generateSpectralShiftVectors(circSizeX, circSizeY, 0, 0);

    forwardX = repmat(forwardX, 1, circSizeY);
    backwardX = repmat(backwardX, 1, circSizeY);
    forwardY = repmat(forwardY, circSizeX, 1);
    backwardY = repmat(backwardY, circSizeX, 1);
end
