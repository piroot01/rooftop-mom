function values = evaluateRemainder(distance, wavenumber)
% evaluateRemainder: Returns the remainder exp(-ikR) / R - S(R) of the Green's function after the extraction of
%   S(R) = 1 / R - k^2 R / 2 + k^4 R^3 / 24 (generateExtractedMoments). It tends to -ik as R tends to zero.
%
%   Inputs:
%       distance: Array of positive distances R.
%       wavenumber: Wavenumber k.
%
%   Output:
%       values: Complex array of the size of distance.

    % expm1 avoids the cancellation in exp(-ikR) - 1 for small kR.
    values = expm1(-1i * wavenumber * distance) ./ distance + wavenumber^2 / 2 * distance ...
        - wavenumber^4 / 24 * distance.^3;
end
