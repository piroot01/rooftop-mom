function [forwardX, backwardX, forwardY, backwardY] = generateSpectralShiftOperators(circSizeX, circSizeY)
% generateSpectralShiftOperators: Compute phase multipliers corresponding to spatial shifts in
%   the Fourier domain. The shift used here is exp(i*pi*u/sizeX)
%   and similarly for Y. The backward shift is the complex
%   conjugate of the forward shift.

    [u, v] = ndgrid(0:(circSizeX - 1), 0:(circSizeY - 1));

    forwardX = exp((2 * 1i * pi * u) / circSizeX);
    forwardY = exp((2 * 1i * pi * v) / circSizeY);

    backwardX = conj(forwardX);
    backwardY = conj(forwardY);
end
