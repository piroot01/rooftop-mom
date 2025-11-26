classdef UniformRectPrimitiveUniformRooftops < rooftopMoM.basisFunctions.UniformRectPrimitiveDivergenceConforming
    properties (Access = private, Constant)
        origin_ = [0; 0; 0];
        eX_ = [1; 0; 0];
        eY_ = [0; 1; 0];
    end

    properties (Access = private)
        quadrature_ rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature;
    end

    methods
        function obj = UniformRectPrimitiveUniformRooftops(primitive, quadrature)
            obj@rooftopMoM.basisFunctions.UniformRectPrimitiveDivergenceConforming(primitive);
            obj.quadrature_ = quadrature;
        end
    end

    methods
        function generate(obj)
            generate@rooftopMoM.basisFunctions.UniformRectPrimitiveDivergenceConforming(obj);
            obj.quadrature_.generate();
        end

        function [quadraturePoints, rooftopPoints, rooftopValues] = generateRooftopsDir(obj, dir)
        % generateRooftopsDirection: Generate rooftop basis for given axis direction.

            % Fetch primitive dimensions and grid size.
            [sX, sY] = obj.primitive().faceSize();
            [cX, cY] = obj.primitive().gridSize();

            % Determine count and shifts for the requested direction.
            switch dir
                case 'x'
                    % Number of rooftops along x-directed edges.
                    count = double((cX - 1) * cY);
                    shifts = rooftopMoM.basisFunctions.utils.generateRectGridShifts(cX - 1, cY, sX, sY);
                case 'y'
                    % Number of rooftops along y-directed edges.
                    count = double((cY - 1) * cX);
                    shifts = rooftopMoM.basisFunctions.utils.generateRectGridShifts(cX, cY - 1, sX, sY);
                otherwise
                    error("Unknown direction '%s'. Use 'x' or 'y'.", direction);
            end

            % Generate a single rooftop's geometry and values.
            [singlePoints, rooftopValues] = obj.generateRooftop_(dir);
            orderSquared = obj.quadrature_.order()^2;
            columnCount = orderSquared * count;

            quadraturePoints = zeros(3, columnCount, 2);
            shiftsDistributed = repelem(shifts, 1, orderSquared);
            quadraturePoints(:, :, 1) = repmat(singlePoints.quadraturePointsPositive, 1, count) + shiftsDistributed;
            quadraturePoints(:, :, 2) = repmat(singlePoints.quadraturePointsNegative, 1, count) + shiftsDistributed;

            rooftopPoints = zeros(3, columnCount, 2);
            rooftopPoints(:, :, 1) = repmat(singlePoints.rooftopPointsPositive, 1, count);
            rooftopPoints(:, :, 2) = repmat(singlePoints.rooftopPointsNegative, 1, count);
        end

        function quadrature = quadrature(obj)
            quadrature = obj.quadrature_;
        end
    end

    methods (Access = private)
        function [rooftopPoints, rooftopValues] = generateRooftop_(obj, direction)
            switch direction
                case 'x'
                    e = obj.eX_;
                    len = obj.primitive().sizeX();
                    ht = obj.primitive().sizeY();
                case 'y'
                    e = obj.eY_;
                    len = obj.primitive().sizeY();
                    ht = obj.primitive().sizeX();
                otherwise
                    error("Unknown direction '%s'. Use 'x' or 'y'.", direction);
            end

            % Define the geometry of the two rooftop wedges by theirs shifts.
            shiftPos = obj.origin_;
            shiftNeg = shiftPos + len * e;
            shifts = cat(2, shiftPos, shiftNeg);

            % Generate the quadrature points.
            e1 = obj.primitive().sizeX() * obj.eX_;
            e2 = obj.primitive().sizeY() * obj.eY_;
            [pointsU, pointsV] = obj.quadrature().points();
            [qpPos, qpNeg] = rooftopMoM.basisFunctions.utils.generateQuadraturePointsUniformRooftop(e1, e2, shifts, pointsU, pointsV);

            % The first rooftop wedge is positive and has an origin in [0; 0; 0].
            rpPos = (1 / len) * (e * (e.' * (qpPos - shiftPos)));

            % The second one is negative and must be shifted by the total length of the rooftop (2 * len) in the e direction.
            rpNeg = (-1 / len) * (e * (e.' * (qpNeg - 2 * shiftNeg)));

            rooftopPoints = struct( ...
                'quadraturePointsPositive', qpPos, ...
                'quadraturePointsNegative', qpNeg, ...
                'rooftopPointsPositive', rpPos, ...
                'rooftopPointsNegative', rpNeg ...
            );

            rooftopValues = struct( ...
                'rooftopDivPositive', 1 / len, ...
                'rooftopDivNegative', (-1) / len, ...
                'length', len, ...
                'height', ht ...
            );
        end
    end
end
