function result = singularity_6(a, b)
    % computes the integral of: l1 * (1 - l1') * R
    sqrt_term = sqrt(a^2 + b^2);

    log_term1 = log((a + sqrt_term) / b);
    log_term2 = log((b + sqrt_term) / a);

    term1 = (sqrt_term * (-90 * a^6 + 213 * a^4 * b^2 - 4 * a^2 * b^4 + 8 * b^6)) / (3780 * a^2);
    term2 = (90 * a^7 - 8 * b^7 + 105 * a^3 * b^4 * log_term1 + 210 * a^6 * b * log_term2) / (3780 * a^2);

    result = term1 + term2;
end