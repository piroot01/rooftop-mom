classdef Mesh < handle
    properties (Access = public)
        % 3 by N matrix of vertex coords
        vertices double;

        % faces will be represented as an array of index lists which represent a vertices of given face
        faces uint32;

        % edges as a pair of vertex indices, the innerEdges will be used for basis function construction
        innerEdges uint32;
        outerEdges uint32;
    end
end