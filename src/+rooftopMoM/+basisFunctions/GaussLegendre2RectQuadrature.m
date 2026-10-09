classdef GaussLegendre2RectQuadrature
% GaussLegendre2RectQuadrature: Tensor Gauss-Legendre rule on the unit square. The B = order^2 points are ordered with
%   the u coordinate running fastest.

    properties (Access = private)
        order_ double;
        pointsU_ double;
        pointsV_ double;
        weights2D_ double;
    end

    methods
        function obj = GaussLegendre2RectQuadrature(order)
            [points, weights] = rooftopMoM.basisFunctions.utils.generateGaussLegendre1D(order, 0, 1);
            [pointsU, pointsV] = ndgrid(points, points);
            [weightsU, weightsV] = ndgrid(weights, weights);

            obj.order_ = order;
            obj.pointsU_ = pointsU(:).';
            obj.pointsV_ = pointsV(:).';
            obj.weights2D_ = weightsU(:) .* weightsV(:);
        end
    end

    methods
        function order = order(obj)
            order = obj.order_;
        end

        function [pointsU, pointsV] = points(obj)
        % points: Returns the 1xB coordinates u and v of the points on the unit square.

            pointsU = obj.pointsU_;
            pointsV = obj.pointsV_;
        end

        function weights = weights2D(obj)
        % weights2D: Returns the Bx1 weights on the unit square.

            weights = obj.weights2D_;
        end
    end
end
