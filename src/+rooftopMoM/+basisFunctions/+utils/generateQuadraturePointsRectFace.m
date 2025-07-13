function points = generateQuadraturePointsRectFace(e1, e2, shift, points2Du, points2Dv)
% generateQuadraturePointsRectFace: Function to compute quadrature points on a rectangular face
%
%   Inputs:
%       e1: 3×1 vector for the first principal edge of the face.
%       e2: 3×1 vector for the second principal edge of the face.
%       shift: 3×1 translation vector for the face origin.
%       points2Du: 1×N array of reference quadrature points along the u‑axis.
%       points2Dv: 1×N array of reference quadrature points along the v‑axis.
%
%   Outputs:
%       points: 3×N array of transformed quadrature points on the face.

    points = e1 .* points2Du + e2 .* points2Dv + shift;
end
