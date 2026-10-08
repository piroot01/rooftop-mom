function [quadraturePoints, rooftopPoints] = convertRooftopPointsToArrays(points)
% convertRooftopPointsToArrays: Convert rooftop point structs into 3D arrays.
%
%   Input:
%       points: struct from UniformRectPrimitiveUniformRooftops.
%
%   Outputs:
%       quadraturePoints: 3xNx2 array of quadrature points.
%       rooftopPoints: 3xNx2 array of rooftop points.

    quadraturePoints = rooftopMoM.core.greensMatrix.utils.stackPositiveAndNegativeFields(points, 'quadraturePointsPositive', 'quadraturePointsNegative');
    rooftopPoints = rooftopMoM.core.greensMatrix.utils.stackPositiveAndNegativeFields(points, 'rooftopPointsPositive', 'rooftopPointsNegative');
end
