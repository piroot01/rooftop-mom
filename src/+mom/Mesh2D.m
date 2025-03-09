classdef Mesh2D < handle
    properties (Access = public)
        % 3 by N matrix of vertex coords
        vertices double;

        % faces will be represented as an array of index lists which represent a vertices of given face
        faces uint32;
    end

    methods
        function vertices = getVertices(obj)
            vertices = obj.vertices;
        end

        function faces = getFaces(obj)
            faces = obj.faces;
        end
    end
end