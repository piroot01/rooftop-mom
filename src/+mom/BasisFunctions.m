classdef BasisFunctions < handle
    properties (Access = private)
        % store the rooftop basis functions as a pair of face indices
        mesh_;

        faces_ uint32;
        funcIDs_ uint32;

        quadPoints_ double;
        rooftopPoints_ double;
        rooftopPointsDiv_ double;

        rooftopLengths_ double;
        rooftopLengthsIdx_ uint32;
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

            function faces = getDoubleFaces(edges)
                sortedEdges = sortrows(edges, [1, 2]);
                diffEdges = diff(sortedEdges(:, 1:2));

                doubleEdgesIndices = find(all(diffEdges == 0, 2));

                face1 = sortedEdges(doubleEdgesIndices, 3);
                face2 = sortedEdges(doubleEdgesIndices + 1, 3);

                faces = [min(face1, face2), max(face1, face2)];

                % TODO: think of a better basis function finding algorithm
                % faces = sortrows(faces, [1, 2]);
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
            obj.faces_ = getDoubleFaces(generateEdges(mesh.getFaces())).';
            obj.funcIDs_ = generateIDs(mesh, obj.faces_);
        end
    end

    methods
        function generatePoints(obj, quadrature)
            % input:
            %   mesh :: mom.Mesh2D
            %   quadrature :: mom.Quadrature

            basisFuncCount = obj.getCount();
            order2 = quadrature.getOrder()^2;

            pointCount = order2 * basisFuncCount;

            obj.quadPoints_ = zeros(3, pointCount, 2);
            obj.rooftopPoints_ = zeros(3, pointCount, 2);
            obj.rooftopPointsDiv_ = zeros(1, basisFuncCount, 2);

            rooftopLengthsMap = containers.Map('KeyType', 'char', 'ValueType', 'uint32');
            obj.rooftopLengths_ = [];
            obj.rooftopLengthsIdx_ = zeros(1, basisFuncCount);

            meshFaces = obj.mesh_.getFaces();
            meshVertices = obj.mesh_.getVertices();

            for n = 1:basisFuncCount
                bIdx = obj.faces_(:, n);
                vIdx = meshFaces(:, bIdx);
                f = meshVertices(:, vIdx(:));
                q = quadrature.get(f);

                [pP, pN, l1, l2] = obj.computeRooftop(f, q);

                slice = (order2 * (n - 1) + 1):(order2 * n);

                obj.quadPoints_(:, slice, 1) = q(:, 1:4);
                obj.quadPoints_(:, slice, 2) = q(:, 5:8);
                obj.rooftopPoints_(:, slice, 1) = pP;
                obj.rooftopPoints_(:, slice, 2) = pN;
                obj.rooftopPointsDiv_(:, n, 1) = 1 / l1;
                obj.rooftopPointsDiv_(:, n, 2) = -1 / l1;

                key = sprintf('%.8f_%.8f', l1, l2);

                if rooftopLengthsMap.isKey(key)
                    obj.rooftopLengthsIdx_(n) = rooftopLengthsMap(key);
                else
                    idx = size(obj.rooftopLengths_, 2) + 1;
                    rooftopLengthsMap(key) = idx;
                    obj.rooftopLengths_(:, end + 1) = [l1; l2];
                    obj.rooftopLengthsIdx_(n) = idx;
                end
            end
        end

        function count = getCount(obj)
            count = size(obj.faces_, 2);
        end

        function faces = getFaces(obj)
            faces = obj.faces_;
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

        function points = getQuadPoints(obj)
            points = obj.quadPoints_;
        end

        function points = getRooftopPoints(obj)
            points = obj.rooftopPoints_;
        end

        function points = getRooftopPointsDiv(obj)
            points = obj.rooftopPointsDiv_;
        end

        function [lengths, indices] = getRooftopLengths(obj)
            lengths = obj.rooftopLengths_;
            indices = obj.rooftopLengthsIdx_;
        end

        function [a, b] = getFaceSize(obj, ID)
            if isa(obj.mesh_, 'mom.Mesh2DCluster')
                [a, b] = obj.mesh_.getMeshSize(ID);
            else
                [a, b] = obj.mesh_.getMeshSize();
            end
        end
    end

    methods (Access = private)
        function [rooftopP, rooftopN, l1, l2] = computeRooftop(~, basisFunctionMesh, quadraturePoints)
            % input:
            %   basisFunctionMesh: 3x8 matrix (8 vertices - technically 6 is enough)
            %   polarity: 0 = positive, 1 = negative
            %   quadraturePoints: 3xN matrix (N is the order of the qudrature squared)

            v1p = basisFunctionMesh(:, 1);
            v3p = basisFunctionMesh(:, 3);
            v1n = basisFunctionMesh(:, 5);
            vn = basisFunctionMesh(:, 7);

            d1 = v1n - v1p;
            l1 = sqrt(sum(d1 .^ 2));
            e = d1 / l1;
            g = v3p - v1p;
            d2 = g - d1;
            l2 = sqrt(sum(d2 .^ 2));

            rooftopP = (1 / l1) * (e * (e.' * (quadraturePoints(:, 1:4) - v1p)));
            rooftopN = (-1 / l1) * (e * (e.' * (quadraturePoints(:, 5:8) - vn)));
        end
    end
end