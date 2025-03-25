classdef Mesh2DCluster < handle
    properties (Access = private)
        meshes_ mom.Mesh2D = mom.Mesh2D.empty;

        lastFaceIndex_ uint32 = 0;
        % starting and ending indexes of each mesh for fast face finding
        % meshRanges uint32 = [];
    end

    methods
        function mesh = getMesh(obj, index)
            mesh = obj.meshes_(index);
        end

        function loadMesh(obj, mesh)
            meshVertexCount = size(mesh.vertices, 2);

            mesh.faces = mesh.faces + obj.lastFaceIndex_;
            obj.meshes_ = [obj.meshes_, mesh];
            obj.lastFaceIndex_ = obj.lastFaceIndex_ + meshVertexCount;
        end

        function vertices = getVertices(obj)
            verticesCell = arrayfun(@(m) m.getVertices(), obj.meshes_, 'UniformOutput', false);
            vertices = [verticesCell{:}];
        end

        function faces = getFaces(obj)
            faceCells = arrayfun(@(m) m.getFaces(), obj.meshes_, 'UniformOutput', false);
            faces = [faceCells{:}];
        end

        function count = getMeshCount(obj)
            count = numel(obj.meshes_);
        end

        %{
        function ID = getMeshID(obj, face)
            ID = find(face >= obj.meshRanges(:, 1) & face <= obj.meshRanges(:, 2), 1, 'first');

            if isempty(ID)
                error('The face was not found!');
            end
        end
        %}

        function [a, b] = getMeshSize(obj, ID)
            mesh = obj.meshes_(ID);
            sampleFace = mesh.faces(:, 1);
            vertices = obj.getVertices();

            v1 = vertices(:, sampleFace(1));
            v2 = vertices(:, sampleFace(2));
            v4 = vertices(:, sampleFace(4));

            a = norm(v2 - v1);
            b = norm(v4 - v1);
        end
    end
end
