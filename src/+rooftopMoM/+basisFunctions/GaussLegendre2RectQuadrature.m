classdef GaussLegendre2RectQuadrature < handle
    properties (Access = private)
        order_ uint32;

        points2Du_ double;
        points2Dv_ double;
        weights2D_ double;
        weights4D_ double;
    end

    methods
        function obj = GaussLegendre2RectQuadrature(order)
            obj.order_ = order;
        end
    end

    methods
        function generate(obj)
            [points1D, weights1D] = rooftopMoM.basisFunctions.utils.generateGaussLegendre1D(double(obj.order_), 0, 1);

            [p2Du, p2Dv] = ndgrid(points1D, points1D);
            obj.points2Du_ = p2Du(:).';
            obj.points2Dv_ = p2Dv(:).';

            [w2D1, w2D2] = ndgrid(weights1D, weights1D);
            obj.weights2D_ = w2D1(:) .* w2D2(:);
            obj.weights4D_ = obj.weights2D_ * obj.weights2D_.';
        end

        function [points2Du, points2Dv] = points(obj)
            points2Du = obj.points2Du_;
            points2Dv = obj.points2Dv_;
        end

        function weights = weights2D(obj)
            weights = obj.weights2D_;
        end

        function weights = weights4D(obj)
            weights = obj.weights4D_;
        end

        function order = order(obj)
            order = obj.order_;
        end
    end
end
