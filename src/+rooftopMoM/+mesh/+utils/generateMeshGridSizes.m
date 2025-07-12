function sizes = generateMeshGridSizes(mesh)
% generateMeshGridSizes: Compute maximal grid sizes for each planar mesh in a cluster.
%
%   Input:
%       mesh: Instance of Mesh or MeshCluster class.
%
%   Output:
%       sizes: Cell array of [c1, c2] pairs, one per mesh in the cluster, where c1 and c2 are grid counts along the two principal axes.

    meshCount = mesh.count();

    [c1s, c2s] = arrayfun( ...
        @(k) rooftopMoM.mesh.utils.generateRectGridSize( ...
            mesh.get(k).vertices(), mesh.get(k).faces() ...
            ), ...
        1:meshCount ...
    );

    % Convert two row vectors into a cell.
    sizes = num2cell([c1s(:), c2s(:)], 2);
end
