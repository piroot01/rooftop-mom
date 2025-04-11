classdef BasisFunctions < handle
    properties (Access = private)
        mesh_;
        isMesh2DCluster_ logical;

        facesX_;
        facesY_;

        funcIDs_ uint32;
    end

    methods
        function obj = BasisFunctions(mesh)
            obj.mesh_ = mesh;

            if isa(mesh, 'mom.Mesh2DCluster')
                obj.isMesh2DCluster_ = true;
                mesh2DCount = mesh.getMeshCount();

                for n = 1:mesh2DCount
                    [facesX, facesY] = obj.generateBasisFunctionFaceIndices(mesh.getMesh(n));
                    obj.facesX_{n} = facesX;
                    obj.facesY_{n} = facesY;
                end
            else
                obj.isMesh2DCluster_ = false;
                [obj.facesX_, obj.facesY_] = obj.generateBasisFunctionFaceIndices(obj.mesh_);
            end
        end
    end

    methods
        function count = getCount(obj)
            if obj.isMesh2DCluster_
                count = size([obj.facesX_{:}], 2) + size([obj.facesY_{:}], 2);
            else
                count = size(obj.facesX_, 2) + size(obj.facesY_, 2);
            end
        end

        function faces = getFaces(obj)
            if obj.isMesh2DCluster_
                faces = obj.getFacesAsCell();
                offsets = cumsum([0, cellfun(@(x) max([x(:); 0]), faces)]);
                faces = cellfun(@(x, offset) x + offset, faces, num2cell(offsets(1:(end - 1))), 'UniformOutput', false);

                faces = [faces{:}];
            else
                faces = [obj.facesX_, obj.facesY_];
            end
        end

        function faces = getFacesAsCell(obj)
            if obj.isMesh2DCluster_
                mesh2DCount = obj.mesh_.getMeshCount();
                facesX = obj.facesX_;
                facesY = obj.facesY_;

                faces = arrayfun(@(n) [facesX{n}, facesY{n}], 1:mesh2DCount, 'UniformOutput', false);
            else
                faces{1} = [obj.facesX_, obj.facesY_];
            end
        end

        function faces = getFacesXAsCell(obj)
            if obj.isMesh2DCluster_
                mesh2DCount = obj.mesh_.getMeshCount();
                facesX = obj.facesX_;

                faces = arrayfun(@(n) [facesX{n}], 1:mesh2DCount, 'UniformOutput', false);
            else
                faces{1} = obj.facesX_;
            end
        end

        function faces = getFacesYAsCell(obj)
            if obj.isMesh2DCluster_
                mesh2DCount = obj.mesh_.getMeshCount();
                facesY = obj.facesY_;

                faces = arrayfun(@(n) [facesY{n}], 1:mesh2DCount, 'UniformOutput', false);
            else
                faces{1} = obj.facesY_;
            end
        end

        function faces = getFacesInDirection(obj, direction)
            if direction == 'x'
                faces = obj.getFacesXAsCell();
            elseif direction == 'y'
                faces = obj.getFacesYAsCell();
            else
                error("Unknown direction: %s\n", direction);
            end
        end

        function IDs = getIDs(obj)
            faces = obj.getFacesAsCell();
            meshCount = numel(faces);

            BasisFunctionCounts = cellfun(@(x) size(x, 2), faces);

            IDs = repelem(1:meshCount, BasisFunctionCounts);
        end

        function count = getIDCount(obj)
            %{
            if obj.isMesh2DCluster_
                count = obj.mesh_.getMeshCount();
            else
                count = 1;
            end
            %}

            % fancier version
            count = obj.isMesh2DCluster_ * obj.mesh_.getMeshCount() + ~obj.isMesh2DCluster_;
        end
    end

    methods (Access = private)
        function [facesX, facesY] = generateBasisFunctionFaceIndices(~, mesh)
            function edges = generateEdges(faces)
                faceCount = size(faces, 2);

                edges = zeros(4 * faceCount, 3, 'uint32');

                faceIndices = (1:faceCount)';

                v1 = faces(1, :);
                v2 = faces(2, :);
                v3 = faces(3, :);
                v4 = faces(4, :);

                edges(1:4:end, :) = [min(v1, v2).', max(v1, v2).', faceIndices];
                edges(2:4:end, :) = [min(v2, v3).', max(v2, v3).', faceIndices];
                edges(3:4:end, :) = [min(v3, v4).', max(v3, v4).', faceIndices];
                edges(4:4:end, :) = [min(v4, v1).', max(v4, v1).', faceIndices];
            end

            function [facesX, facesY] = getDoubleFaces(edges)
                sortedEdges = sortrows(edges, [1, 2]);
                diffEdges = diff(sortedEdges(:, 1:2));

                doubleEdgesIndices = find(all(diffEdges == 0, 2));

                face1 = sortedEdges(doubleEdgesIndices, 3);
                face2 = sortedEdges(doubleEdgesIndices + 1, 3);

                faces = [min(face1, face2), max(face1, face2)];

                faceDiff = faces(:, 2) - faces(:, 1);
                adjacentIdx = (faceDiff == 1);

                facesX = faces(adjacentIdx, :);
                facesY = faces(~adjacentIdx, :);
            end

            [facesX, facesY] = getDoubleFaces(generateEdges(mesh.getFaces()));
            facesX = uint32(facesX.');
            facesY = uint32(facesY.');
        end
    end
end
