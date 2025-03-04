classdef Mesh2DCluster < handle
    properties (Access = public)
        meshes mom.Mesh2D = mom.Mesh2D.empty;

        % starting and ending indexes of each mesh for fast face finding
        meshRanges uint32 = [];
    end

    methods
        function mesh = getMesh(obj, index)
            mesh = obj.meshes(index);
        end

        function loadMesh(obj, mesh)
            function vertexIndex = getLastIndex(mesh)
                vertexIndex = uint32(size(mesh(end).vertices, 2));
            end

            currMeshLastIndex = getLastIndex(mesh);
            if isempty(obj.meshes)
                obj.meshRanges = [obj.meshRanges; 1, currMeshLastIndex];
            else
                prevMeshLastIndex = getLastIndex(obj.meshes(end));
                obj.meshRanges = [obj.meshRanges; prevMeshLastIndex + 1, prevMeshLastIndex + currMeshLastIndex];
            end

            mesh.faces = mesh.faces + currMeshLastIndex;
            obj.meshes = [obj.meshes, mesh];
        end

        function count = getMeshCount(obj)
            count = numel(obj.meshes);
        end
    end
end