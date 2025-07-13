classdef MeshInterface < handle
    methods (Abstract)
        mesh = get(obj, index);

        vertices = vertices(obj);

        faces = faces(obj);

        count = count(obj);
    end
end
