function sizes = generateGridSizes(mesh)
% generateGridSizes: Computes the maximal rectangular grid size to fit the plannar mesh in.
%
%   Input:
%       mesh: Instance of a class Mesh or MeshCluster.
%   Output:
%       sizes: Cell of pairs, each pair correcponds to one mesh grid sizes.

    meshCount = mesh.count();

    if meshCount == 1
        vertices = mesh.vertices();

        [coeff, score] = pca(vertices.')

        arrayfun()
        min(norm(vertices))
    end
end
