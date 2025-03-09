function result = singularity_2(a, b)
    % computes the integral of: l1 * l1' * R
    sqrt_term = sqrt(a^2 + b^2);
    log_term1 = log((a + sqrt_term) / b);
    log_term2 = log((b + sqrt_term) / a);

    term1 = sqrt_term * (-36 * a^6 + 165 * a^4 * b^2 - 122 * a^2 * b^4 - 8 * b^6);
    term2 = (36 * a^7 + 126 * a^2 * b^5 + 8 * b^7 + ...
            210 * a^3 * b^4 * log_term1 + 105 * a^6 * b * log_term2);

    result = term1 / (3780 * a^2) + term2 / (3780 * a^2);
end