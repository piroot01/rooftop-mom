function [length, height] = generateRooftopParametersUniformRooftop(support)
% generateRooftopParametersUniformRooftop: Computes the metrics of a rooftop on specified support.
%
%   Inputs:
%       support: A 3×8 array of the eight vertices defining the two adjacent parallelograms:
%                [v1p, v2p, v3p, v4p,   % Omega_i^+ corners
%                 v1n, v2n, v3n, v4n]   % Omega_i^– corners
%
%   Outputs:
%       length: Half of the total length of the longer side of the rooftop.
%       height: Size of the shorter side.

    v1p = support(:, 1);
    v3p = support(:, 3);
    v1n = support(:, 5);

    e1 = v1n - v1p;
    length = norm(e1);
    height = norm((v3p - v1p) - e1);
end
