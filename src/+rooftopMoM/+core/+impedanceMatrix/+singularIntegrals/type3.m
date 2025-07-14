function result = type3(a, b)
% type3: Computes the integral of 1 / R.
%
% Inputs:
%   a: Half the length of the longer sides of the rooftop.
%   b: Length of the shorter side of the rooftop.
%
% Output:
%   result: Analytical value of the singular integral.

    sqrt_term = sqrt(a.^2 + b.^2);
    log_term1 = log((a + sqrt_term) ./ b);
    log_term2 = log((b + sqrt_term) ./ a);

    numerator = 2 * (a.^3 - a.^2 .* sqrt_term + b.^3 - b.^2 .* sqrt_term ...
                     + 3 * a .* b.^2 .* log_term1 + 3 * a.^2 .* b .* log_term2);

    result = numerator ./ (3 * a.^2 .* b.^2);
end

