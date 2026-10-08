function result = type1(a, b)
% type1: Computes the integral of l1 * l1' / R.
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

    numerator = 12 * a.^5 - 12 * a.^4 .* sqrt_term ...
                + a.^2 .* b.^2 .* (30 * b - 29 * sqrt_term) ...
                + 2 * b.^4 .* (b - sqrt_term) ...
                + 15 * a.^3 .* b .* (4 * b .* log_term1 + 3 * a .* log_term2);

    result = numerator ./ (90 * a.^4 .* b.^2); 
end
