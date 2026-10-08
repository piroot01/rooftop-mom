function value = type2(a, b)
% type2: Returns the integral of u * (1 - u') / R over a face paired with itself, the static part of a rising and a
%   falling rooftop half on the same face. All coordinates run over the unit square and
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

    numerator = 18 * a.^5 - 2 * b.^5 ...
        - 18 * a.^4 .* diagonal ...
        - a.^2 .* b.^2 .* diagonal ...
        + 2 * b.^4 .* diagonal ...
        + 30 * a.^3 .* b.^2 .* logA ...
        + 45 * a.^4 .* b .* logB;

    value = numerator ./ (90 * a.^4 .* b.^2);
end
