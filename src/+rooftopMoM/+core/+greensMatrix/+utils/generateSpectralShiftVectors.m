function [forwardX, backwardX, forwardY, backwardY] = generateSpectralShiftVectors(sizeX, sizeY, parityX, parityY)
% generateSpectralShiftVectors: Computes the phase multipliers of unit spatial shifts in the Fourier domain
%   as vectors ready for implicit expansion. The forward shift at the frequency index u is
%   exp(1i * pi * (2 * u + parity) / size), the backward shift is its complex conjugate.
%
%   Inputs:
%       sizeX: Positive integer specifying the spectrum size in the x direction.
%       sizeY: Positive integer specifying the spectrum size in the y direction.
%       parityX: 0 for the integer frequencies, 1 for the frequencies shifted by half a bin in the x direction.
%       parityY: 0 for the integer frequencies, 1 for the frequencies shifted by half a bin in the y direction.
%
%   Outputs:
%       forwardX: sizeX x 1 vector of the forward shift in the x direction.
%       backwardX: sizeX x 1 vector of the backward shift in the x direction.
%       forwardY: 1 x sizeY vector of the forward shift in the y direction.
%       backwardY: 1 x sizeY vector of the backward shift in the y direction.

    forwardX = exp(1i * pi * (2 * (0:(sizeX - 1)).' + parityX) / sizeX);
    forwardY = exp(1i * pi * (2 * (0:(sizeY - 1)) + parityY) / sizeY);

    backwardX = conj(forwardX);
    backwardY = conj(forwardY);
end
