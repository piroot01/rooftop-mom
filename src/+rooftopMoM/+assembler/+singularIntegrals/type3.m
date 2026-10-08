function value = type3(a, b)
% type3: Returns the integral of 1 / R over a face paired with itself, the static part of the scalar potential on
%   the same face. All coordinates run over the unit square and R = sqrt((a * (u - u'))^2 + (b * (v - v'))^2).
%
%   Inputs:
%       a: Face size in the x direction (u).
%       b: Face size in the y direction (v).
%
%   Output:
%       value: Value of the integral.

    diagonal = sqrt(a.^2 + b.^2);
    logA = log((a + diagonal) ./ b);
    logB = log((b + diagonal) ./ a);

    numerator = 2 * (a.^3 - a.^2 .* diagonal + b.^3 - b.^2 .* diagonal ...
        + 3 * a .* b.^2 .* logA + 3 * a.^2 .* b .* logB);

    value = numerator ./ (3 * a.^2 .* b.^2);
end
