classdef MeshInterface < handle
    methods (Abstract)
        mesh = get(obj, index);

        vertices = vertices(obj);

        faces = faces(obj);

        count = count(obj);

        memory = memory(obj);
    end
end
