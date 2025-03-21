classdef RooftopBasisFunctions < handle
    properties (Access = private)
        mesh_;

        basisFunctionIndices_ uint32;

        quadraturePoints_ double;
        rooftopPoints_ double;
        rooftopDivPoints_ double;

        uniqueRooftopSizes_ double;
        uniqueRooftopSizeIndices_ uint32;
    end

    methods
        function obj = RooftopBasisFunctions(mesh, basisFunctionIndices)
            obj.mesh_ = mesh;
            obj.basisFunctionIndices_ = basisFunctionIndices;
        end
    end

    methods
        function generatePoints(obj, quadrature)
            basisFunctionIndices = obj.basisFunctionIndices_;
            basisFunctionCount = size(basisFunctionIndices, 2);

            orderSquared = quadrature.getOrder()^2;

            pointCount = orderSquared * basisFunctionCount;

            obj.quadraturePoints_ = zeros(3, pointCount, 2);
            obj.rooftopPoints_ = zeros(3, pointCount, 2);
            obj.rooftopDivPoints_ = zeros(1, pointCount, 2);

            meshVertices = obj.mesh_.getVertices();
            meshFaces = obj.mesh_.getFaces();

            allRooftopSizes = zeros(2, basisFunctionCount);

            for n = 1:basisFunctionCount
                singleBasisFunctionIndices = basisFunctionIndices(:, n);
                basisFunctionFaces = meshFaces(:, singleBasisFunctionIndices);
                basisFunctionVertices = meshVertices(:, basisFunctionFaces(:));

                basisFunctionQuadraturePoints = quadrature.get(basisFunctionVertices);

                [rooftopPointsPositive, rooftopPointsNegative, l1, l2] = obj.computeRooftopValues(basisFunctionVertices, basisFunctionQuadraturePoints);
                allRooftopSizes(:, n) = [l1; l2];

                slice = ((n - 1) * orderSquared + 1):(n * orderSquared);
                obj.quadraturePoints_(:, slice, 1) = basisFunctionQuadraturePoints(:, 1:4); 
                obj.quadraturePoints_(:, slice, 2) = basisFunctionQuadraturePoints(:, 5:8); 
                obj.rooftopPoints_(:, slice, 1) = rooftopPointsPositive;
                obj.rooftopPoints_(:, slice, 2) = rooftopPointsNegative;
                obj.rooftopDivPoints_(:, n, 1) = 1 / l1;
                obj.rooftopDivPoints_(:, n, 2) = (-1) / l1;
            end

            % NOTE: the round is limiting the resolution of the overall algorithm precision, but 10 places are fine :D
            [obj.uniqueRooftopSizes_, ~, obj.uniqueRooftopSizeIndices_] = unique(round(allRooftopSizes.', 10), 'rows', 'stable');
        end

        function points = getQuadraturePoints(obj)
            points = obj.quadraturePoints_;
        end

        function points = getRooftopPoints(obj)
            points = obj.rooftopPoints_;
        end

        function points = getRooftopDivPoints(obj)
            points = obj.rooftopDivPoints_;
        end

        function [sizes, indices] = getUniqueRooftopSizes(obj)
            sizes = obj.uniqueRooftopSizes_;
            indices = obj.uniqueRooftopSizeIndices_;
        end
    end

    methods (Access = private)
        function [rooftopP, rooftopN, l1, l2] = computeRooftopValues(~, basisFunctionMesh, quadraturePoints)
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
