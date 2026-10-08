function points = generateSingleQuadraturePointsUniformRectPrimitive(quadrature, primitive)
% generateSingleQuadraturePointsUniformRectPrimitive: Returns quadrature points for all faces of a uniform rectangular primitive.
%
%   Inputs:
%       quadrature: GaussLegendre2RectQuadrature instance (provides points() and order()).
%       primitive: UniformRectPrimitive instance (provides gridSize() and faceSize()).
%
%   Output:
%       points: 3 x N array of 3D points, where
%           N = quadrature.order()^2

    % Unit basis vectors (3x1).
    eX = [1; 0; 0];
    eY = [0; 1; 0];

    % Single face (cell) dimensions in x and y.
    [sX, sY] = primitive.faceSize();

    % Edge vectors for a single rectangular face of the primitive.
    e1 = sX * eX;
    e2 = sY * eY;

    % Quadrature points on a reference unit face (2D).
    [points2Du, points2Dv] = quadrature.points();

    % Generate quadrature points for a single origin face (3 x blockSize).
    points = rooftopMoM.basisFunctions.utils.generateQuadraturePointsRectFace( ...
        e1, e2, [0; 0; 0], points2Du, points2Dv);
end
