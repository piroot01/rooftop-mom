classdef BasisFunctions < handle
    properties (Access = private)
        % store the rooftop basis functions as a pair of face indices
        faces_ uint32;

        quadPoints_ double;
        rooftopPoints_ double;
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

            obj.faces_ = getDoubleFaces(generateEdges(mesh.getFaces())).';
        end
    end

    methods
        function generatePoints(obj, mesh, quadrature)
            % input:
            %   mesh :: mom.Mesh2D
            %   quadrature :: mom.Quadrature

            basisFuncCount = obj.getCount();
            order2 = quadrature.getOrder()^2;

            pointCount = order2 * basisFuncCount;

            obj.quadPoints_ = zeros(3, pointCount, 2);
            obj.rooftopPoints_ = zeros(3, pointCount, 2);

            for n = 1:basisFuncCount
                bIdx = obj.faces_(:, n);
                vIdx = mesh.faces(:, bIdx);
                f = mesh.vertices(:, vIdx(:));
                q = quadrature.get(f);
                [pP, pN] = obj.computeRooftop(f, q);
                slice = (order2 * (n - 1) + 1):(order2 * n);
                obj.quadPoints_(:, slice, 1) = q(:, 1:4);
                obj.quadPoints_(:, slice, 2) = q(:, 5:8);
                obj.rooftopPoints_(:, slice, 1) = pP;
                obj.rooftopPoints_(:, slice, 2) = pN;
            end
        end

        function count = getCount(obj)
            count = size(obj.faces_, 2);
        end

        function faces = getFaces(obj)
            faces = obj.faces_;
        end

        function points = getQuadPoints(obj)
            points = obj.quadPoints_;
        end

        function points = getRooftopPoints(obj)
            points = obj.rooftopPoints_;
        end
    end

    methods (Access = private)
        function [rooftopP, rooftopN] = computeRooftop(~, basisFunctionMesh, quadraturePoints)
            % input:
            %   basisFunctionMesh: 3x8 matrix (8 vertices - technically 6 is enough)
            %   polarity: 0 = positive, 1 = negative
            %   quadraturePoints: 3xN matrix (N is the order of the qudrature squared)

            v1p = basisFunctionMesh(:, 1);
            v1n = basisFunctionMesh(:, 5);
            v2n = basisFunctionMesh(:, 6);

            % Compute the difference vector and its norm.
            d1 = v1n - v1p;
            d2 = v2n - v1n;
            l1 = sqrt(sum(d1 .^ 2));
            l2 = sqrt(sum(d2 .^ 2));
            e1 = d1 / l1;
            e2 = d2 / l2;

            % Compute the rooftop using a vectorized projection.
            rooftopP = (1 / l1) * (e1 * (e1.' * (quadraturePoints(:, 1:4) - v1p)));
            rooftopN = (-1 / l2) * (e2 * (e2.' * (quadraturePoints(:, 5:8) - v2n)));
        end
    end
end