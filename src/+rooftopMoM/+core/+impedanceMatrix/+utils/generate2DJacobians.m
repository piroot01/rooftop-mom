function jacobians = generate2DJacobians(mesh, IDMask)
% generate2DJacobians: Returns a 1D array of Jacobians for all uniform rooftops.
%
%   Inputs:
%       mesh: Structure containing vertices and faces for multiple meshes.
%       IDMask: Array mapping each rooftop to a mesh index.
%
%   Output:
%       jacobians: 1xN array, where each element is the Jacobian for 2D integration using single rooftop.

    % Compute one area per rooftop in a vector.
    areas = cellfun( ...
       @(vs, fs) rooftopMoM.mesh.utils.generateRectFaceArea(vs, fs), ...
       mesh.vertices, mesh.faces );

    % Then just look up each mask‐entry in that vector.
    jacobians = areas(IDMask);
end
