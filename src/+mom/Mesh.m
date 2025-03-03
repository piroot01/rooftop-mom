classdef Mesh
    properties (Access = private)
        vertices Vec3;

        % faces will be represented as an array of index lists which represent a vertices of given face
        faces int32;

        % edges as a pair of vertex indices, the innerEdges will be used for basis function construction
        innerEdges int32;
        outerEdges int32;
    end
end