classdef RooftopBasisFunctions < handle
    properties (Access = private, Constant)
        ex_ = [1; 0; 0];
        ey_ = [0; 1; 0];
    end

    properties (Access = private)
        quadrature_ mom.Quadrature;
        gridSize_ uint32;
    end

    methods
        function obj = RooftopBasisFunctions(quadrature, gridSize)
            obj.quadrature_ = quadrature;
            obj.gridSize_ = gridSize;
        end
    end

    methods
        function gridSize = getGridSize(obj)
            gridSize = obj.gridSize_;
        end

        function count = getCount(obj)
            count = prod(obj.gridSize_);
        end

        function points = getRooftopPointsInDirection(obj, direction)
            basisFunctionCount = prod(obj.gridSize_);
            quadraturePoints = obj.quadrature_.getQuadraturePointsRectangle(obj.ex_, obj.ey_);

            if direction == 'x'
                [rooftopPointsPlusXDir, rooftopPointsMinusXDir] = obj.computeRooftopPoints(quadraturePoints, 'x');
                points = cat(3, repmat(rooftopPointsPlusXDir, 1, basisFunctionCount), repmat(rooftopPointsMinusXDir, 1, basisFunctionCount));
            elseif direction == 'y'
                [rooftopPointsPlusYDir, rooftopPointsMinusYDir] = obj.computeRooftopPoints(quadraturePoints, 'y');
                points = cat(3, repmat(rooftopPointsPlusYDir, 1, basisFunctionCount), repmat(rooftopPointsMinusYDir, 1, basisFunctionCount));
            else
                error("Unknown direction: %s\n", direction);
            end
        end

        function points = getQuadraturePointsInDirection(obj, direction)
            quadraturePoints = obj.quadrature_.getQuadraturePointsRectangle(obj.ex_, obj.ey_);

            orderSquared = obj.quadrature_.getOrder()^2;

            countX = obj.gridSize_(1) + 1;
            countY = obj.gridSize_(2) + 1;

            [mGrid, nGrid] = ndgrid(1:countX, 1:countY);

            shift = cat(2, mGrid(:) - 1, nGrid(:) - 1, zeros(countX * countY, 1));
            quadraturePointsReshaped = reshape(quadraturePoints, 3, orderSquared, 1, 1);
            shiftReshaped = reshape(shift.', 3, 1, countX, countY);

            allQuadraturePoints = double(shiftReshaped) + quadraturePointsReshaped;

            if direction == 'x'
                indicesXPlus = 1:(obj.gridSize_(1));
                indicesYPlus = 1:(obj.gridSize_(2));
                indicesXMinus = 2:(obj.gridSize_(1) + 1);
                indicesYMinus = 1:(obj.gridSize_(2));

                points = cat(3, ...
                    reshape(allQuadraturePoints(:, :, indicesXPlus, indicesYPlus), 3, []), ...
                    reshape(allQuadraturePoints(:, :, indicesXMinus, indicesYMinus), 3, []));
            elseif direction == 'y'
                indicesXPlus = 1:(obj.gridSize_(1));
                indicesYPlus = 1:(obj.gridSize_(2));
                indicesXMinus = 1:(obj.gridSize_(1));
                indicesYMinus = 2:(obj.gridSize_(2) + 1);

                points = cat(3, ...
                    reshape(allQuadraturePoints(:, :, indicesXPlus, indicesYPlus), 3, []), ...
                    reshape(allQuadraturePoints(:, :, indicesXMinus, indicesYMinus), 3, []));
            elseif direction == 'a'
                points = allQuadraturePoints;
            else
                error("Unknown direction: %s\n", direction);
            end
        end

        function points = getQuadraturePoints(obj, padded)
            if padded
                points = obj.getQuadraturePointsInDirection('a');
                return
            end

            quadraturePoints = obj.quadrature_.getQuadraturePointsRectangle(obj.ex_, obj.ey_);
            orderSquared = size(quadraturePoints, 2);

            Nx = obj.gridSize_(1);
            Ny = obj.gridSize_(2);

            numX = (Nx - 1) * Ny;
            numY = Nx * (Ny - 1);
            numCols = orderSquared * (numX + numY);

            points = zeros(3, numCols, 2, 'like', quadraturePoints);

            function slice = buildSlice(iList, jList)
                [I,J] = ndgrid(iList, jList);
                shifts = [I(:)'; J(:)'; zeros(1,numel(I))];    
                tmp = bsxfun(@plus, reshape(quadraturePoints, 3, orderSquared, 1), reshape(double(shifts), 3, 1, []));
                slice = reshape(tmp, 3, []);
            end

            XP = buildSlice(0:(Nx - 2), 0:(Ny - 1));
            YP = buildSlice(0:(Nx - 1), 0:(Ny - 2));
            points(:,:,1) = [XP, YP];

            XM = buildSlice(1:(Nx - 1), 0:(Ny - 1));
            YM = buildSlice(0:(Nx - 1), 1:(Ny - 1));
            points(:,:,2) = [XM, YM];
        end

        %{
        function points = getQuadraturePoints(obj, padded)
            if padded == true
                points = obj.quadraturePoints_;
            else
                indicesXPlusXDir = 1:(obj.gridSize_(1) - 1);
                indicesXPlusYDir = 1:obj.gridSize_(2);

                indicesXMinusXDir = 2:obj.gridSize_(1);
                indicesXMinusYDir = indicesXPlusYDir;

                indicesYPlusXDir = 1:obj.gridSize_(1);
                indicesYPlusYDir = 1:(obj.gridSize_(2) - 1);

                indicesYMinusXDir = indicesYPlusXDir;
                indicesYMinusYDir = 2:obj.gridSize_(2);

                XP = reshape(obj.quadraturePoints_(:, :, indicesXPlusXDir, indicesXPlusYDir), 3, []);
                XM = reshape(obj.quadraturePoints_(:, :, indicesXMinusXDir, indicesXMinusYDir), 3, []);
                YP = reshape(obj.quadraturePoints_(:, :, indicesYPlusXDir, indicesYPlusYDir), 3, []);
                YM = reshape(obj.quadraturePoints_(:, :, indicesYMinusXDir, indicesYMinusYDir), 3, []);

                points = cat(3, cat(2, XP, YP), cat(2, XM, YM));
            end
        end
        %}

        function points = getRooftopPoints(obj)
            quadraturePoints = obj.quadrature_.getQuadraturePointsRectangle(obj.ex_, obj.ey_);

            [rooftopPointsPlusXDir, rooftopPointsMinusXDir] = obj.computeRooftopPoints(quadraturePoints, 'x');
            [rooftopPointsPlusYDir, rooftopPointsMinusYDir] = obj.computeRooftopPoints(quadraturePoints, 'y');

            gridSizeXDir = obj.gridSize_(1);
            gridSizeYDir = obj.gridSize_(2);

            basisFunctionCountXDir = (gridSizeXDir - 1) * gridSizeYDir;
            basisFunctionCountYDir = gridSizeXDir * (gridSizeYDir - 1);

            rooftopPoitsXDir = cat(3, repmat(rooftopPointsPlusXDir, 1, basisFunctionCountXDir), repmat(rooftopPointsMinusXDir, 1, basisFunctionCountXDir));
            rooftopPoitsYDir = cat(3, repmat(rooftopPointsPlusYDir, 1, basisFunctionCountYDir), repmat(rooftopPointsMinusYDir, 1, basisFunctionCountYDir));

            points = cat(2, rooftopPoitsXDir, rooftopPoitsYDir);
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
