function result = type2(a, b)
    % computes the integral of: l1 * (1 - l1') / R
    assert(isequal(size(a), size(b)), 'Inputs a and b must have the same dimensions!');

    sqrt_term = sqrt(a.^2 + b.^2);
    log_term1 = log((a + sqrt_term) ./ b);
    log_term2 = log((b + sqrt_term) ./ a);

    numerator = 18 * a.^5 - 2 * b.^5 ...
                - 18 * a.^4 .* sqrt_term ...
                - a.^2 .* b.^2 .* sqrt_term ...
                + 2 * b.^4 .* sqrt_term ...
                + 30 * a.^3 .* b.^2 .* log_term1 ...
                + 45 * a.^4 .* b .* log_term2;

    result = numerator ./ (90 * a.^4 .* b.^2);
end