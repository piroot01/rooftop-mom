function offsets = generateMeshOffsets(mesh)
% generateMeshOffsets: Retrieve the offsets based on the face counts of each mesh plane.
%
%   Input:
%       mesh: Instance of Mesh or MeshCluster class.
%
%   Output:
%       offsets: 1xN array of face offsets

    meshCount = mesh.count();

    if meshCount == 1
        offsets = 0;
    else
        faces = mesh.faces();
        offsets = uint32([0; arrayfun(@(n) size(faces{n}, 2), 2:meshCount, 'UniformOutput', true)]);
    end
end
