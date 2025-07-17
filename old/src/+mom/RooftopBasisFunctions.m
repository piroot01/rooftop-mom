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
            obj.rooftopDivPoints_ = zeros(1, basisFunctionCount, 2);

            meshVertices = obj.mesh_.getVertices();
            meshFaces = obj.mesh_.getFaces();

            allRooftopSizes = zeros(2, basisFunctionCount);

            for n = 1:basisFunctionCount
                singleBasisFunctionIndices = basisFunctionIndices(:, n);
                basisFunctionFaces = meshFaces(:, singleBasisFunctionIndices);
                basisFunctionVertices = meshVertices(:, basisFunctionFaces(:));

                basisFunctionQuadraturePoints = quadrature.get(basisFunctionVertices);

                [rooftopPointsPositive, rooftopPointsNegative, l1, l2] = obj.computeRooftopValues(basisFunctionVertices, basisFunctionQuadraturePoints, orderSquared);
                allRooftopSizes(:, n) = [l1; l2];

                slice = ((n - 1) * orderSquared + 1):(n * orderSquared);
                obj.quadraturePoints_(:, slice, 1) = basisFunctionQuadraturePoints(:, 1:orderSquared); 
                obj.quadraturePoints_(:, slice, 2) = basisFunctionQuadraturePoints(:, (orderSquared + 1):end); 
                obj.rooftopPoints_(:, slice, 1) = rooftopPointsPositive;
                obj.rooftopPoints_(:, slice, 2) = rooftopPointsNegative;
                obj.rooftopDivPoints_(:, n, 1) = 1 / l1;
                obj.rooftopDivPoints_(:, n, 2) = (-1) / l1;
            end

            % NOTE: the round is limiting the resolution of the overall algorithm precision, but 10 places are fine :D
            [obj.uniqueRooftopSizes_, ~, obj.uniqueRooftopSizeIndices_] = unique(round(allRooftopSizes.', 10), 'rows', 'stable');
        end

        function points = getQuadraturePoints(obj, ~)
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

        function sizes = getRooftopSizes(obj)
            sizes = obj.uniqueRooftopSizes_(obj.uniqueRooftopSizeIndices_, :);
        end

        function sizes = getUniqueUniqueRooftopSizes(obj, IDs)
            indices = [1, find(diff(IDs)) + 1];
            sizes = obj.uniqueRooftopSizes_(obj.uniqueRooftopSizeIndices_(indices), :);
        end

        % temporary - just for testing
        function points = getWedgePoints(~, quadraturePoints, alpha, polarity)
            if alpha == 'x'
                e = [1; 0; 0];
            else
                e = [0; 1; 0];
            end

            if polarity == 'p'
                points = e * (e.' * quadraturePoints);
            else
                points = -1 * e * (e.' * (quadraturePoints - e));
            end
        end
    end

    methods (Access = private)
        function [rooftopP, rooftopN, l1, l2] = computeRooftopValues(~, basisFunctionMesh, quadraturePoints, blockSize)
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

            rooftopP = (1 / l1) * (e * (e.' * (quadraturePoints(:, 1:blockSize) - v1p)));
            rooftopN = (-1 / l1) * (e * (e.' * (quadraturePoints(:, (blockSize + 1):end) - vn)));
        end

    end
end
