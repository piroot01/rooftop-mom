function [pointsPositive, pointsNegative] = generateQuadraturePointsUniformRooftop(e1, e2, shifts, points2Du, points2Dv)
% generateQuadraturePointsUniformRooftop: Assemble quadrature points for both sides of a uniform rooftop element.
%
%   Inputs:
%       e1: 3×1 vector for the first principal edge of the face.
%       e2: 3×1 vector for the second principal edge of the face.
%       shifts: 3×2 array of translation vectors [shift_pos, shift_neg].
%       points2Du: 1×N array of reference quadrature points along the u‑axis.
%       points2Dv: 1×N array of reference quadrature points along the v‑axis.
%
%   Outputs:
%       pointsPositive: 3×N array of transformed quadrature points for positive faces.
%       pointsNegative: 3×N array of transformed quadrature points for negative faces.

    % Positive side quadrature.
    pointsPositive = rooftopMoM.basisFunctions.utils.generateQuadraturePointsRectFace( ...
        e1, e2, shifts(:,1), points2Du, points2Dv);

    % Negative side quadrature.
    pointsNegative = rooftopMoM.basisFunctions.utils.generateQuadraturePointsRectFace( ...
        e1, e2, shifts(:,2), points2Du, points2Dv);
end
