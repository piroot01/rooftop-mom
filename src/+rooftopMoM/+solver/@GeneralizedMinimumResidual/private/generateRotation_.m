function [cosine, sine] = generateRotation_(a, b)
% generateRotation_: Returns the complex Givens rotation [conj(cosine), conj(sine); -sine, cosine] that maps [a; b]
%   onto [norm([a, b]); 0].
%
%   Inputs:
%       a, b: Complex scalars.
%
%   Outputs:
%       cosine, sine: Complex scalars, abs(cosine)^2 + abs(sine)^2 = 1.

    radius = norm([a, b]);
    cosine = 1;
    sine = 0;

    if radius > 0
        cosine = a / radius;
        sine = b / radius;
    end
end
