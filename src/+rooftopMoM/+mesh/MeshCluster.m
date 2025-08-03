classdef MeshCluster < rooftopMoM.mesh.MeshInterface
    properties (Access = private)
        meshes_ cell = {};
    end

    methods
        function mesh = get(obj, index)
            mesh = obj.meshes_{index};
        end

        function insert(obj, mesh)
            obj.meshes_{end + 1} = mesh;
        end

        function remove(obj, index)
            obj.meshes_{index} = [];
        end

        function vertices = vertices(obj)
            vertices = cellfun(@(mesh) [mesh.vertices{1}], obj.meshes_, 'UniformOutput', false);
        end

        function faces = faces(obj)
            faces = cellfun(@(mesh) [mesh.faces{1}], obj.meshes_, 'UniformOutput', false);
        end

        function count = count(obj)
            count = numel(obj.meshes_);
        end

        function memory = memory(obj)
            bytes = sum(cellfun(@(mesh) mesh.memory().bytes(), obj.meshes_));
            memory = rooftopMoM.Memory(bytes, 'B');
        end
    end
end
