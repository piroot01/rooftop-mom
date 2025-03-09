function result = singularity_1(a, b)
    % computes the integral of: l1 * l1' / R
    assert(isequal(size(a), size(b)), 'Inputs a and b must have the same dimensions!');

    sqrt_term = sqrt(a.^2 + b.^2);
    log_term1 = log((a + sqrt_term) ./ b);
    log_term2 = log((b + sqrt_term) ./ a);

    term1 = sqrt_term .* (-12*a.^4 - 29*(a.^2 .* b.^2) - 2*b.^4);
    term2 = 12*a.^5 + 30*(a.^2 .* b.^3) + 2*b.^5 + ...
            60*(a.^3 .* b.^2 .* log_term1) + 45*(a.^4 .* b .* log_term2);

    result = (term1 + term2) ./ (90 * a.^2);
end