classdef Mesh < handle
    properties (Access = private)
        % 3×N array where each column is the [x; y; z] coordinates of one vertex.
        vertices_ double;

        % F×N array where each column lists the indices of the F vertices that make up one face.
        faces_ uint32;
    end

    methods
        function obj = Mesh(vertices, faces)
            obj.vertices_ = vertices;
            obj.faces_ = faces;
        end
    end

    methods
        function obj = get(obj, ~)
        end

        function vertices = vertices(obj)
            vertices = obj.vertices_;
        end

        function faces = faces(obj)
            faces = obj.faces_;
        end

        function setVertices(obj, vertices)
            obj.vertices_ = vertices;
        end

        function setFaces(obj, faces)
            obj.faces_ = faces;
        end

        function count = count(~)
            count = uint32(1);
        end
    end
end
