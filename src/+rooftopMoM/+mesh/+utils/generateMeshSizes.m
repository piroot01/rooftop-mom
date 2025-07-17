function sizes = generateMeshSizes(mesh)
% generateMeshSizes: Compute maximal sizes for each planar mesh in a cluster.
%
%   Input:
%       mesh: Instance of Mesh or MeshCluster class.
%
%   Output:
%       sizes: Cell array of [s1, s2] pairs, one per mesh in the cluster, where s1 and s2 are sizes along the two principal axes.

    meshCount = mesh.count();

    [s1s, s2s] = arrayfun( ...
        @(k) rooftopMoM.mesh.utils.generateRectSize( ...
            mesh.get(k).vertices{1}, mesh.get(k).faces{1} ...
            ), ...
        1:meshCount ...
    );

    % Convert two row vectors into a cell.
    sizes = num2cell([s1s(:), s2s(:)], 2);
end
