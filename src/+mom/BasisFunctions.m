classdef BasisFunctions < handle
    properties (Access = private)
        % store the rooftop basis functions as a pair of face indices
        faces_ uint32;

        % store the rooftops (positive and negative)
        pointsPositive_ double;
        pointsNegative_ double;
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
            end

            obj.faces_ = getDoubleFaces(generateEdges(mesh.getFaces())).';
        end
    end

    methods
        function generatePoints(obj, mesh, quadraturePoints)
            % input:
            %   mesh :: Mesh2D
            %   quadraturePoints :: QuadPoints

            basisFuncCount = obj.getCount();
            obj.pointsPositive_ = zeros(3, basisFuncCount);
            order2 = double(quadraturePoints.getOrder())^2;

            obj.pointsPositive_ = zeros(3, order2 * basisFuncCount);
            obj.pointsNegative_ = zeros(3, order2 * basisFuncCount);

            for n = 1:basisFuncCount
                bIdx = obj.faces_(:, n);
                vIdx = mesh.faces(:, bIdx);
                f = mesh.vertices(:, vIdx(:));
                q = quadraturePoints.getTransformedBasisFunc(f);
                [pP, pN] = obj.computeRooftop(f, q);
                slice = (order2 * (n - 1) + 1):(order2 * n);
                obj.pointsPositive_(:, slice) = pP;
                obj.pointsNegative_(:, slice) = pN;
            end
        end

        function count = getCount(obj)
            count = size(obj.faces_, 2);
        end

        function [pP, pN] = getPoints(obj)
            pP = obj.pointsPositive_;
            pN = obj.pointsNegative_;
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