function sizes = generateMeshSizes(mesh)
% generateMeshSizes: Compute maximal sizes for each planar mesh in a cluster.
%
%   Input:
%       mesh: Instance of Mesh or MeshCluster class.
%
%   Output:
%       sizes: Cell array of [s1, s2] pairs, one per mesh in the cluster, where s1 and s2 are sizes along the two principal axes.

    meshCount = mesh.count();
    sizes = cell(1, meshCount);

    for n = 1:meshCount
        vertices = mesh.get(n).vertices{1};
        faces = mesh.get(n).faces{1};

        [u, v] = rooftopMoM.mesh.utils.generateRectLocalCoords(vertices, faces);

        [s1, s2] = rooftopMoM.mesh.utils.generateRectSize(vertices, u, v);

        sizes{n} = [s1, s2];
    end
end
