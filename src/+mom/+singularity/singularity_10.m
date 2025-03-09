function result = singularity_10(a, b)
    % computes the integral of: R
    sqrt_term = sqrt(a^2 + b^2);
    log_term1 = log((a + sqrt_term) / b);
    log_term2 = log((b + sqrt_term) / a);

    term1 = (1/60) * sqrt_term * (-4 * a^4 + 12 * a^2 * b^2 - 4 * b^4);
    term2 = (1/60) * (4 * a^5 + 4 * b^5 + 10 * a * b^4 * log_term1 + 10 * a^4 * b * log_term2);

    result = term1 + term2;
end