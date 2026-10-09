function values = integrateMonomial(power, m, n, a, b)
% integrateMonomial: Returns the mixed antiderivative G(a, b) of the near-field note in docs/nearfield, the integral
%   of u^m v^n (u^2 + v^2)^(power / 2) over [0, a] x [0, b], in its closed form (Closed forms, Appendix A): with
%   R = sqrt(a^2 + b^2), integer polynomials P0 to P3 and an integer D,
%       G = (P0 + P1 R + P2 log((R + b) / a) + P3 log((R + a) / b)) / D.
%   G vanishes on both axes, so the integral over a rectangle of the first quadrant is the alternating sum of G at its
%   corners.
%
%   Inputs:
%       power: -1, 1 or 3.
%       m: 0, 1, 2 or 3.
%       n: 0 or 1.
%       a, b: Nonnegative arrays of one size.
%
%   Output:
%       values: Array of the size of a.

    R = sqrt(a.^2 + b.^2);

    % log((R + b) / a) and log((R + a) / b), as asinh to keep the accuracy where the quotient is close to one.
    logB = asinh(b ./ a);
    logA = asinh(a ./ b);

    switch sprintf('%d,%d,%d', power, m, n)
        case '-1,0,0'
            values = a .* logB + b .* logA;
        case '-1,0,1'
            values = (-a.^2 + a .* R + b.^2 .* logA) / 2;
        case '-1,1,0'
            values = (-b.^2 + b .* R + a.^2 .* logB) / 2;
        case '-1,1,1'
            values = (-a.^3 - b.^3 + (a.^2 + b.^2) .* R) / 3;
        case '-1,2,0'
            values = (a .* b .* R + 2 * a.^3 .* logB - b.^3 .* logA) / 6;
        case '-1,2,1'
            values = (-2 * a.^4 + a .* (2 * a.^2 + b.^2) .* R - b.^4 .* logA) / 8;
        case '-1,3,0'
            values = (2 * b.^4 + b .* (a.^2 - 2 * b.^2) .* R + 3 * a.^4 .* logB) / 12;
        case '-1,3,1'
            values = (-3 * a.^5 + 2 * b.^5 + (3 * a.^4 + a.^2 .* b.^2 - 2 * b.^4) .* R) / 15;
        case '1,0,0'
            values = (2 * a .* b .* R + a.^3 .* logB + b.^3 .* logA) / 6;
        case '1,0,1'
            values = (-2 * a.^4 + a .* (2 * a.^2 + 5 * b.^2) .* R + 3 * b.^4 .* logA) / 24;
        case '1,1,0'
            values = (-2 * b.^4 + b .* (5 * a.^2 + 2 * b.^2) .* R + 3 * a.^4 .* logB) / 24;
        case '1,1,1'
            values = (-a.^5 - b.^5 + (a.^4 + 2 * a.^2 .* b.^2 + b.^4) .* R) / 15;
        case '1,2,0'
            values = (a .* b .* (6 * a.^2 + b.^2) .* R + 4 * a.^5 .* logB - b.^5 .* logA) / 40;
        case '1,2,1'
            values = (-8 * a.^6 + a .* (8 * a.^4 + 14 * a.^2 .* b.^2 + 3 * b.^4) .* R - 3 * b.^6 .* logA) / 144;
        case '1,3,0'
            values = (4 * b.^6 + b .* (21 * a.^4 + 2 * a.^2 .* b.^2 - 4 * b.^4) .* R + 15 * a.^6 .* logB) / 180;
        case '1,3,1'
            values = (-5 * a.^7 + 2 * b.^7 + (5 * a.^6 + 8 * a.^4 .* b.^2 + a.^2 .* b.^4 - 2 * b.^6) .* R) / 105;
        case '3,0,0'
            values = (7 * a .* b .* (a.^2 + b.^2) .* R + 3 * a.^5 .* logB + 3 * b.^5 .* logA) / 40;
        case '3,0,1'
            values = (-8 * a.^6 + a .* (8 * a.^4 + 26 * a.^2 .* b.^2 + 33 * b.^4) .* R + 15 * b.^6 .* logA) / 240;
        case '3,1,0'
            values = (-8 * b.^6 + b .* (33 * a.^4 + 26 * a.^2 .* b.^2 + 8 * b.^4) .* R + 15 * a.^6 .* logB) / 240;
        case '3,1,1'
            values = (-a.^7 - b.^7 + (a.^6 + 3 * a.^4 .* b.^2 + 3 * a.^2 .* b.^4 + b.^6) .* R) / 35;
        case '3,2,0'
            values = (a .* b .* (38 * a.^4 + 26 * a.^2 .* b.^2 + 3 * b.^4) .* R + 18 * a.^7 .* logB ...
                - 3 * b.^7 .* logA) / 336;
        case '3,2,1'
            values = (-48 * a.^8 + a .* (48 * a.^6 + 136 * a.^4 .* b.^2 + 118 * a.^2 .* b.^4 + 15 * b.^6) .* R ...
                - 15 * b.^8 .* logA) / 1920;
        case '3,3,0'
            values = (16 * b.^8 + b .* (215 * a.^6 + 134 * a.^4 .* b.^2 + 8 * a.^2 .* b.^4 - 16 * b.^6) .* R ...
                + 105 * a.^8 .* logB) / 2240;
        case '3,3,1'
            values = (-7 * a.^9 + 2 * b.^9 ...
                + (7 * a.^8 + 19 * a.^6 .* b.^2 + 15 * a.^4 .* b.^4 + a.^2 .* b.^6 - 2 * b.^8) .* R) / 315;
    end

    % Every term tends to zero on the axes, where a logarithm term evaluates to 0 * Inf or 0 / 0.
    values(a == 0 | b == 0) = 0;
end
