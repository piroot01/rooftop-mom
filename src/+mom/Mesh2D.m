classdef Mesh2D < handle
    properties (Access = public)
        % 3 by N matrix of vertex coords
        vertices double;

        % faces will be represented as an array of index lists which represent a vertices of given face
        faces uint32;

        % mesh plane description
        % origin must be a 3D column vector
        origin double;
        transformation double;
    end

    methods
        function faces = getFaces(obj)
            faces = obj.faces;
        end
    end
end