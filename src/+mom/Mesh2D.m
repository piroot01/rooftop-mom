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

        function [a, b] = getMeshSize(obj)
            sampleFace = obj.faces(:, 1);

            v1 = obj.vertices(:, sampleFace(1));
            v2 = obj.vertices(:, sampleFace(2));
            v4 = obj.vertices(:, sampleFace(4));

            a = norm(v2 - v1);
            b = norm(v4 - v1);
        end
    end
end
