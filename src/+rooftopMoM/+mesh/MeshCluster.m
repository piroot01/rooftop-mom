classdef MeshCluster < handle
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
            verticesCell = cellfun(@(mesh) mesh.vertices(), obj.meshes_, 'UniformOutput', false);
            vertices = [verticesCell{:}];
        end

        function faces = faces(obj)
            meshCount = obj.count();

            facesCell = cell(1, meshCount);
            offset = uint32(0);

            for n = 1:meshCount
                vertices = obj.get(n).vertices();
                faces = obj.get(n).faces();

                % Offset this mesh’s face indices so they point into the global list.
                facesCell{n} = faces + offset;

                % Increase offset by how many vertices this mesh contributed.
                offset = offset + uint32(size(vertices, 2));
            end

            faces = [facesCell{:}];
        end

        function count = count(obj)
            count = numel(obj.meshes_);
        end
    end
end
