classdef QuadPoints
    properties
        order_ uint32;
        points2Du_ double;
        points2Dv_ double;
        weights4D_ double;
    end

    methods
        function obj = QuadPoints(order)
            function [points, weights] = generateQuadPoints1D(order)
                [points, weights] = mom.lgwt(double(order), 0, 1);
            end

            obj.order_ = order;

            [p1D, w1D] = generateQuadPoints1D(order);

            [p2Du, p2Dv] = ndgrid(p1D, p1D);
            obj.points2Du_ = p2Du(:).';
            obj.points2Dv_ = p2Dv(:).';

            [w2D1, w2D2] = ndgrid(w1D, w1D);
            weights2D = w2D1(:) .* w2D2(:);
            obj.weights4D_ = weights2D * weights2D.';
        end
    end

    methods
        function points = getTransformed(obj, faceVertices)
            % we expect the faceVertices to be 3 by (4N) matrix, where N is the number of faces
            faceCount = size(faceVertices, 2) / 4;

            v1 = faceVertices(:, 1:4:end);
            v2 = faceVertices(:, 2:4:end);
            v4 = faceVertices(:, 4:4:end);

            uSide= v2 - v1;
            vSide = v4 - v1;

            uCount = numel(obj.points2Du_);
            points = zeros(3, faceCount * uCount);

            for n = 1:faceCount
                slice = ((n - 1) * uCount + 1):(n * uCount);
                facePoints = uSide(:, n) .* obj.points2Du_ + vSide(:, n) .* obj.points2Dv_ + v1(:, n);
                points(:, slice) = facePoints;
            end
        end

        function points = getTransformedBasisFunc(obj, basisFunctionMesh)
            % each basis function consists of two faces
            faceCount = 2;

            v1 = basisFunctionMesh(:, 1:4:end);
            v2 = basisFunctionMesh(:, 2:4:end);
            v4 = basisFunctionMesh(:, 4:4:end);

            uSide= v2 - v1;
            vSide = v4 - v1;

            uCount = numel(obj.points2Du_);
            points = zeros(3, faceCount * uCount);

            facePointsP = uSide(:, 1) .* obj.points2Du_ + vSide(:, 1) .* obj.points2Dv_ + v1(:, 1);
            facePointsN = uSide(:, 2) .* obj.points2Du_ + vSide(:, 2) .* obj.points2Dv_ + v1(:, 2);

            points(:, 1:uCount) = facePointsP;
            points(:, (uCount + 1):end) = facePointsN;
        end

        function order = getOrder(obj)
            order = obj.order_;
        end

        function weights = getWeights4D(obj)
            weights = obj.weights4D_;
        end
    end
end