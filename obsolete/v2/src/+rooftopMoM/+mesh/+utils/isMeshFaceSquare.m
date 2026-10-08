function verdict = isMeshFaceSquare(mesh)
% isMeshFaceSquare: Determine whether given mesh has square mesh.
%
%   Input:
%       mesh: Instance of Mesh or MeshCluster class.
%
%   Output:
%       verdict: Logical value indicating whether the specified face is square.

    meshCount = mesh.count();

    verdict = arrayfun( ...
        @(k) rooftopMoM.mesh.utils.isRectFaceSquare( ...
            mesh.get(k).vertices{1}, mesh.get(k).faces{1} ...
            ), ...
        1:meshCount ...
    );
end
