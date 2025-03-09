function result = singularity_9(a, b)
    % computes the integral of: 1 / R
    sqrt_term = sqrt(a^2 + b^2);
    log_term1 = log((a + sqrt_term) / b);
    log_term2 = log((b + sqrt_term) / a);

    term1 = -(2/3) * (a^2 + b^2) * sqrt_term;
    term2 = (2/3) * (a^3 + b^3 + 3 * a * b^2 * log_term1 + 3 * a^2 * b * log_term2);

    result = term1 + term2;
end