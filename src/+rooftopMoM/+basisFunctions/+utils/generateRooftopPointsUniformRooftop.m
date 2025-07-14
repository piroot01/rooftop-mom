function [pointsPositive, pointsNegative] = generateRooftopPointsUniformRooftop(support, quadraturePointsPositive, quadraturePointsNegative)
% generateRooftopPointsUniformRooftop: Returns the rooftop points for a uniform rooftop.
%
%   Inputs:
%       support: A 3×8 array of the eight vertices defining the two adjacent parallelograms:
%                [v1p, v2p, v3p, v4p,   % Omega_i^+ corners
%                 v1n, v2n, v3n, v4n]   % Omega_i^– corners
%       quadraturePointsPositive: 3×N array of quadrature points on the positive side of the rooftop.
%       quadraturePointsNegative: 3×N array of quadrature points on the negative side of the rooftop.
%
%   Outputs:
%       pointsPositive: 3×N array of projected quadrature points on the positive side.
%       pointsNegative: 3×N array of projected quadrature points on the negative side.

    % Extract the key vertices:  
    %   v1p defines the origin and length of the positive side,
    %   v1n defines the opposite endpoint of that side,
    %   v3n defines the origin of the negative side.
    v1p = support(:, 1);
    v1n = support(:, 5);
    v3n = support(:, 7);

    % Compute the rooftop edge vector’s length and its unit direction u.
    e1 = v1n - v1p;
    l1 = norm(e1);
    u = e1 / l1;

    pointsPositive = (1 / l1) * (u * (u.' * (quadraturePointsPositive - v1p)));
    pointsNegative = (-1 / l1) * (u * (u.' * (quadraturePointsNegative - v3n)));
end
