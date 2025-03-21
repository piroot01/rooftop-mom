classdef BasisFunctions < handle
    properties (Access = private)
        mesh_;

        facesX_ uint32;
        facesY_ uint32;

        funcIDs_ uint32;
    end

    methods
        function obj = BasisFunctions(mesh)
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

            function IDs = generateIDs(mesh, faces)
                if isa(mesh, 'mom.Mesh2DCluster')
                    basisFuncsCount = size(faces, 2);
                    IDs = nan(1, basisFuncsCount);

                    % NOTE: can be further optimised cause the basisFuncs are sorted
                    for n = 1:basisFuncsCount
                        IDs(n) = mesh.getMeshID(faces(1, n));
                    end
                else
                    IDs = ones(1, size(faces, 2));
                end
            end

            obj.mesh_ = mesh;
            [facesX, facesY] = getDoubleFaces(generateEdges(mesh.getFaces()));
            [obj.facesX_, obj.facesY_] = deal(facesX.', facesY.');
            obj.funcIDs_ = generateIDs(mesh, obj.getFaces());
        end
    end

    methods
        function count = getCount(obj)
            count = size(obj.getFaces(), 2);
        end

        function faces = getFaces(obj)
            faces = [obj.facesX_, obj.facesY_];
        end

        function faces = getFacesX(obj)
            faces = obj.facesX_;
        end

        function faces = getFacesY(obj)
            faces = obj.facesY_;
        end

        function IDs = getFuncIDs(obj)
            IDs = obj.funcIDs_;
        end

        function count = getIDCount(obj)
            if isa(obj.mesh_, 'mom.Mesh2DCluster')
                count = obj.mesh_.getMeshCount();
            else
                count = 1;
            end
        end

        function [a, b] = getFaceSize(obj, ID)
            if isa(obj.mesh_, 'mom.Mesh2DCluster')
                [a, b] = obj.mesh_.getMeshSize(ID);
            else
                [a, b] = obj.mesh_.getMeshSize();
            end
        end
    end
end
