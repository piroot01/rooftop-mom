function mask = generateIDMask(basisFunction)
% generateIDMask: Returns the ID mask, each of the basis function index-pair is given an integer based on the mesh to which the function belongs to.
%
%   Input:
%       basisFunction: Any class that derives from UniformRectDivergenceConforming
%
%   Output:
%       mask: 1xM matrix, where M = F_n * N, F_n is the face count of n-th mesh and N is the mesh count

    indices = basisFunction.indicesCell(false);
    meshCount = basisFunction.mesh().count();

    counts = cellfun(@(xs) size(xs, 2), indices);

    mask = repelem(1:meshCount, counts);
end
