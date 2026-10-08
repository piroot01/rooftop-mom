function jacobians = generate2DJacobians(primitive, count)
% generate2DJacobians: Returns a 1D array of Jacobians for all uniform rooftops.
%
%   Inputs:
%       primitive: Structure describing a rectangular antenna plate with no cavities.
%       count: Number of basis functions.
%
%   Output:
%       jacobians: 1xN array, where each element is the Jacobian for 2D integration using single rooftop.

    % Compute one area for single face of the radiator.
    area = primitive.sizeX() * primitive.sizeY();

    % Then just look up each mask‐entry in that vector.
    jacobians = repmat(area, count, 1);
end
