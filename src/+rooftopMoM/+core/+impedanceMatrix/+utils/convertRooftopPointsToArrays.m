function [quadraturePoints, rooftopPoints] = convertRooftopPointsToArrays(points)
% convertRooftopPointsToArrays: Convert rooftop point structs into 3D arrays.
%
%   Input:
%       points: 1×M cell array of structs from UniformRooftops.
%
%   Outputs:
%       quadraturePoints: 3xNx2 array of quadrature points.
%       rooftopPoints: 3xNx2 array of rooftop points.

    quadraturePoints = rooftopMoM.core.impedanceMatrix.utils.stackPositiveAndNegativeFields(points, 'quadraturePointsPositive', 'quadraturePointsNegative');
    rooftopPoints = rooftopMoM.core.impedanceMatrix.utils.stackPositiveAndNegativeFields(points, 'rooftopPointsPositive', 'rooftopPointsNegative');
end
