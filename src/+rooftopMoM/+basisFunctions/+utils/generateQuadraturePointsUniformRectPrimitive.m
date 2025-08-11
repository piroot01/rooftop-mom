function points = generateQuadraturePointsUniformRectPrimitive(quadrature, primitive)
% generateQuadraturePointsUniformRectPrimitive: Returns quadrature points for all faces of a uniform rectangular primitive.
%
%   Inputs:
%       quadrature: GaussLegendre2RectQuadrature instance (provides points() and order()).
%       primitive: UniformRectPrimitive instance (provides gridSize() and faceSize()).
%
%   Output:
%       points: 3 x N array of 3D points, where
%           N = (quadrature.order()^2) * (number of faces).

    % Unit basis vectors (3x1).
    eX = [1; 0; 0];
    eY = [0; 1; 0];

    % Grid dimensions: number of cells in x and y.
    [cX, cY] = primitive.gridSize();

    % Single face (cell) dimensions in x and y.
    [sX, sY] = primitive.faceSize();

    % Edge vectors for a single rectangular face of the primitive.
    e1 = sX * eX;
    e2 = sY * eY;

    % Quadrature points on a reference unit face (2D).
    [points2Du, points2Dv] = quadrature.points();

    % Generate quadrature points for a single origin face (3 x blockSize).
    singlePoints = rooftopMoM.basisFunctions.utils.generateQuadraturePointsRectFace( ...
        e1, e2, [0; 0; 0], points2Du, points2Dv);

    % Grid shifts: positions of each face origin in the full array (3 x faceCount).
    shifts = rooftopMoM.basisFunctions.utils.generateRectGridShifts(cX, cY, sX, sY);

    % Number of quadrature points per face.
    blockSize = quadrature.order()^2;

    % singlePoints: 3 x blockSize
    % shifts: 3 x faceCount, where faceCount = cX * cY
    % We want: all combinations of singlePoints + each shift -> 3 x (blockSize * faceCell)

    % Expand to 3 x blockSize x faceCount then reshape to 3 x (blockSize*faceCount)
    faceCount = size(shifts, 2);

    % reshape to 3 x blockSize x 1 and 3 x 1 x nCells then add
    singlePointsReshaped = reshape(singlePoints, 3, blockSize, 1); % 3 x blockSize x 1
    shiftsReshaped = reshape(shifts, 3, 1, faceCount);  % 3 x 1 x faceCount
    points  = singlePointsReshaped + shiftsReshaped; % implicit expansion: 3 x blockSize x faceCount

    points = reshape(points, 3, []); % 3 x (blockSize * faceCount)
end
