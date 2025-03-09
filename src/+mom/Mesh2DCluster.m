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
            function index = getLastIndex(mesh)
                index = uint32(size(mesh.vertices, 2));
            end

            function faceCount = getFaceCount(mesh)
                faceCount = uint32(size(mesh.faces, 2));
            end

            currFacecount = getFaceCount(mesh);
            if isempty(obj.meshes)
                prevMeshLastIndex = 0;
                obj.meshRanges = uint32([1, currFacecount]);
            else
                prevMeshLastIndex = getLastIndex(obj.meshes(end));
                prevFaceSum = obj.meshRanges(end, 2);
                obj.meshRanges = [obj.meshRanges; prevFaceSum + 1, prevFaceSum + currFacecount];
            end

            mesh.faces = mesh.faces + prevMeshLastIndex;
            obj.meshes = [obj.meshes, mesh];
        end

        function vertices = getVertices(obj)
            meshCount = obj.getMeshCount();
            verticesCells = cell(1, meshCount);

            for n = 1:meshCount
                verticesCells{n} = obj.meshes(n).getVertices();
            end

            vertices = [verticesCells{:}];
        end

        function faces = getFaces(obj)
            meshCount = obj.getMeshCount();
            faceCells = cell(1, meshCount);

            for n = 1:meshCount
                faceCells{n} = obj.meshes(n).getFaces();
            end

            faces = [faceCells{:}];
        end

        function count = getMeshCount(obj)
            count = numel(obj.meshes);
        end

        function ID = getMeshID(obj, face)
            ID = find(face >= obj.meshRanges(:, 1) & face <= obj.meshRanges(:, 2), 1, 'first');

            if isempty(ID)
                error('The face was not found!');
            end
        end
    end
end