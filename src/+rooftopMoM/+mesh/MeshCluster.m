classdef MeshCluster < handle
    properties (Access = private)
        meshes_ cell = {};

        lastOffset_ uint32 = 0;
    end

    methods
        function mesh = get(obj, index)
            mesh = obj.meshes_{index};
        end

        function insert(obj, mesh)
            mesh.setOffset(obj.lastOffset_);
            obj.meshes_{end + 1} = mesh;
            obj.lastOffset_ = obj.lastOffset_ + uint32(size(mesh.vertices(), 2));
        end

        function remove(obj, index)
            obj.lastOffset_ = obj.lastOffset_ - uint32(size(obj.get(index).vertices(), 2));
            obj.meshes_{index} = [];
        end

        function vertices = vertices(obj)
            verticesCell = cellfun(@(mesh) mesh.vertices(), obj.meshes_, 'UniformOutput', false);
            vertices = [verticesCell{:}];
        end

        function faces = faces(obj)
            facesCell = obj.facesCell();
            faces = [facesCell{:}];
        end

        function faces = facesCell(obj)
            faces = arrayfun(@(n) obj.get(n).faces() + obj.get(n).offset(), 1:obj.count(), 'UniformOutput', false);
        end

        function count = count(obj)
            count = numel(obj.meshes_);
        end

        function result = isCluster(~)
            result = true;
        end

        function offset = offset(obj)
            offset = arrayfun(@(n) obj.get(n).offset(), 1:obj.count(), 'UniformOutput', true);
        end
    end
end
