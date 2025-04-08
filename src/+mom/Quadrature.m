classdef Quadrature
    properties (Access = private)
        order_ double;

        points2Du_ double;
        points2Dv_ double;
        weights4D_ double;
        weights2D_ double;
    end

    methods
        function obj = Quadrature(order)
            function [points2Du, points2Dv, weights2D, weights4D] = generate4D(order)
                function [points, weights] = generate1D(order)
                    [points, weights] = mom.lgwt(double(order), 0, 1);
                end

                [p1D, w1D] = generate1D(order);

                [p2Du, p2Dv] = ndgrid(p1D, p1D);
                points2Du = p2Du(:).';
                points2Dv = p2Dv(:).';

                [w2D1, w2D2] = ndgrid(w1D, w1D);
                weights2D = w2D1(:) .* w2D2(:);
                weights4D = weights2D * weights2D.';
            end

            [obj.points2Du_, obj.points2Dv_, obj.weights2D_, obj.weights4D_] = generate4D(order);
            obj.order_ = order;
        end
    end

    methods
        function points = get(obj, basisFunctionMesh)
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

        function weights = getWeights(obj)
            weights = obj.weights4D_;
        end

        function weights = getWeights2D(obj)
            weights = obj.weights2D_;
        end

        function points = getQuadraturePointsRectangle(obj, u, v)
            points = u .* obj.points2Du_ + v .* obj.points2Dv_;
        end
    end
end
