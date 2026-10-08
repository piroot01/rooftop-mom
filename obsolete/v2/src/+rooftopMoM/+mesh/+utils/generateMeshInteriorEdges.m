function indices = generateMeshInteriorEdges(mesh)
% generateMeshInteriorEdges: Compute all interior edges of a mesh (mesh cluster).
%
%   Input:
%       mesh: Instance of Mesh or MeshCluster class.
%
%   Output:
%       indices: Cell of 2xE arrays of face-pairs sharing the interior edge.

    faces = mesh.faces();

    indices = cellfun( ...
        @(xs) rooftopMoM.mesh.utils.generateInteriorEdges(xs), ...
        faces, ...
        'UniformOutput', false ...
    );
end
