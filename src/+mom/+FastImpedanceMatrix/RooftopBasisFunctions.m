classdef RooftopBasisFunctions < handle
    properties (Access = private, Constant)
        ex_ = [1; 0; 0];
        ey_ = [0; 1; 0];
    end

    properties (Access = private)
        quadrature_ mom.Quadrature;
        gridSize_ uint32;

        rooftopPointsXDir_ double;
        rooftopPointsYDir_ double;
    end

    methods
        function obj = RooftopBasisFunctions(quadrature, gridSize)
            obj.quadrature_ = quadrature;
            obj.gridSize_ = gridSize;
        end
    end

    methods
        function generatePoints(obj)
            basisFunctionCount = prod(obj.gridSize_);
            quadraturePoints = obj.quadrature_.getQuadraturePointsRectangle(obj.ex_, obj.ey_);

            [rooftopPointsPlusXDir, rooftopPointsMinusXDir] = obj.computeRooftopPoints(quadraturePoints, 'x');
            [rooftopPointsPlusYDir, rooftopPointsMinusYDir] = obj.computeRooftopPoints(quadraturePoints, 'y');

            obj.rooftopPointsXDir_ = cat(3, repmat(rooftopPointsPlusXDir, 1, basisFunctionCount), repmat(rooftopPointsMinusXDir, 1, basisFunctionCount));
            obj.rooftopPointsYDir_ = cat(3, repmat(rooftopPointsPlusYDir, 1, basisFunctionCount), repmat(rooftopPointsMinusYDir, 1, basisFunctionCount));
        end

        function points = getRooftopPointsInDirection(obj, direction)
            if direction == 'x'
                points = obj.rooftopPointsXDir_;
            elseif direction == 'y'
                points = obj.rooftopPointsYDir_;
            else
                error("Unknown direction: %s\n", direction);
            end
        end
    end

    methods (Access = private)
        function [rooftopPointsPlus, rooftopPointsMinus] = computeRooftopPoints(obj, quadraturePoints, direction)
            if direction == 'x'
                e = obj.ex_;
            elseif direction == 'y'
                e = obj.ey_;
            else
                error("Unknown direction: %s\n", direction);
            end

            rooftopPointsPlus = e * (e.' * quadraturePoints);
            rooftopPointsMinus = (-1) * e * (e.' * (quadraturePoints - e));
        end
    end
end
