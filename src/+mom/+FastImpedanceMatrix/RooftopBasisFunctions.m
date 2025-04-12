classdef RooftopBasisFunctions < handle
    properties (Access = private, Constant)
        ex_ = [1; 0; 0];
        ey_ = [0; 1; 0];
    end

    properties (Access = private)
        quadrature_ mom.Quadrature;
        gridSize_ uint32;

        quadraturePoints_ double;

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

            % generate quadrature points
            orderSquared = obj.quadrature_.getOrder()^2;

            countX = obj.gridSize_(1) + 1;
            countY = obj.gridSize_(2) + 1;

            [mGrid, nGrid] = ndgrid(1:countX, 1:countY);

            shift = cat(2, mGrid(:) - 1, nGrid(:) - 1, zeros(countX * countY, 1));
            quadraturePointsReshaped = reshape(quadraturePoints, 3, orderSquared, 1, 1);
            shiftReshaped = reshape(shift.', 3, 1, countX, countY);

            obj.quadraturePoints_ = double(shiftReshaped) + quadraturePointsReshaped;

            % generate rooftop points
            [rooftopPointsPlusXDir, rooftopPointsMinusXDir] = obj.computeRooftopPoints(quadraturePoints, 'x');
            [rooftopPointsPlusYDir, rooftopPointsMinusYDir] = obj.computeRooftopPoints(quadraturePoints, 'y');

            obj.rooftopPointsXDir_ = cat(3, repmat(rooftopPointsPlusXDir, 1, basisFunctionCount), repmat(rooftopPointsMinusXDir, 1, basisFunctionCount));
            obj.rooftopPointsYDir_ = cat(3, repmat(rooftopPointsPlusYDir, 1, basisFunctionCount), repmat(rooftopPointsMinusYDir, 1, basisFunctionCount));
        end

        function gridSize = getGridSize(obj)
            gridSize = obj.gridSize_;
        end

        function count = getCount(obj)
            count = prod(obj.gridSize_);
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

        function points = getQuadraturePointsInDirection(obj, direction)
            if direction == 'x'
                indicesXPlus = 1:(obj.gridSize_(1));
                indicesYPlus = 1:(obj.gridSize_(2));
                indicesXMinus = 2:(obj.gridSize_(1) + 1);
                indicesYMinus = 1:(obj.gridSize_(2));

                points = cat(3, ...
                    reshape(obj.quadraturePoints_(:, :, indicesXPlus, indicesYPlus), 3, []), ...
                    reshape(obj.quadraturePoints_(:, :, indicesXMinus, indicesYMinus), 3, []));
            elseif direction == 'y'
                indicesXPlus = 1:(obj.gridSize_(1));
                indicesYPlus = 1:(obj.gridSize_(2));
                indicesXMinus = 1:(obj.gridSize_(1));
                indicesYMinus = 2:(obj.gridSize_(2) + 1);

                points = cat(3, ...
                    reshape(obj.quadraturePoints_(:, :, indicesXPlus, indicesYPlus), 3, []), ...
                    reshape(obj.quadraturePoints_(:, :, indicesXMinus, indicesYMinus), 3, []));
            else
                error("Unknown direction: %s\n", direction);
            end
        end

        function points = getQuadraturePoints(obj)
            points = reshape(obj.quadraturePoints_, 3, []);
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
