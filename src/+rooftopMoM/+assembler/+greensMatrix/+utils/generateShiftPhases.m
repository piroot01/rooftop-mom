function [shiftX, shiftY] = generateShiftPhases(sizeX, sizeY, parityX, parityY)
% generateShiftPhases: Returns the phase of a unit spatial shift at the frequencies of a grid,
%   exp(1i * pi * (2 * m + parity) / size), m = 0, ..., size - 1. Parity 1 gives the frequencies halfway between the
%   integer ones, those of the odd subband.
%
%   Inputs:
%       sizeX: Positive integer specifying the grid size in the x direction.
%       sizeY: Positive integer specifying the grid size in the y direction.
%       parityX: 0 or 1, the parity of the frequencies in the x direction.
%       parityY: 0 or 1, the parity of the frequencies in the y direction.
%
%   Outputs:
%       shiftX: sizeX x 1 phases in the x direction.
%       shiftY: 1 x sizeY phases in the y direction.

    shiftX = exp(1i * pi * (2 * (0:(sizeX - 1)).' + parityX) / sizeX);
    shiftY = exp(1i * pi * (2 * (0:(sizeY - 1)) + parityY) / sizeY);
end
