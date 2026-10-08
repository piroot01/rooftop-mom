function value = type1(a, b)
% type1: Returns the integral of u * u' / R over a face paired with itself, the static part of two rising (or two
%   falling) rooftop halves on the same face. All coordinates run over the unit square and
%   R = sqrt((a * (u - u'))^2 + (b * (v - v'))^2).
%
%   Inputs:
%       a: Face size along the rooftop direction (u).
%       b: Face size across the rooftop direction (v).
%
%   Output:
%       value: Value of the integral.

    diagonal = sqrt(a.^2 + b.^2);
    logA = log((a + diagonal) ./ b);
    logB = log((b + diagonal) ./ a);

    numerator = 12 * a.^5 - 12 * a.^4 .* diagonal ...
        + a.^2 .* b.^2 .* (30 * b - 29 * diagonal) ...
        + 2 * b.^4 .* (b - diagonal) ...
        + 15 * a.^3 .* b .* (4 * b .* logA + 3 * a .* logB);

    value = numerator ./ (90 * a.^4 .* b.^2);
end
