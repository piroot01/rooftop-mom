classdef FullRectUniformRooftops < handle
    properties (Access = private, Constant)
        origin_ = [0; 0; 0];
        ex_ = [1; 0; 0];
        ey_ = [0; 1; 0];
    end

    properties (Access = private)
        primitive_ rooftopMoM.mesh.UniformRectPrimitive;

        quadrature_ rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature;
    end

    methods
        function obj = FullRectUniformRooftops(primitive, quadrature)
            obj.primitive_ = primitive;
            obj.quadrature_ = quadrature;
        end
    end

    methods
        function generate(obj)
            obj.quadrature_.generate();
        end

        function rooftopValues = generateRooftops(obj)
            obj.generateRooftop_('y')
        end

        function quadrature = quadrature(obj)
            quadrature = obj.quadrature_;
        end
    end

    methods (Access = private)
        function [rooftopPoints, rooftopValues] = generateRooftop_(obj, direction)
            if direction == 'x'
                e = obj.ex_;
                len = obj.primitive_.sizeX();
                ht = obj.primitive_.sizeY();
            elseif direction == 'y'
                e = obj.ey_;
                len = obj.primitive_.sizeY();
                ht = obj.primitive_.sizeX();
            else
                error("Unknown direction: %s", direction);
            end

            % Define the geometry of the two rooftop wedges by theirs shifts.
            shiftPos = obj.origin_;
            shiftNeg = shiftPos + len * e;
            shifts = cat(2, shiftPos, shiftNeg)

            % Generate the quadrature points.
            e1 = obj.primitive_.sizeX() * obj.ex_;
            e2 = obj.primitive_.sizeY() * obj.ey_;
            [pointsU, pointsV] = obj.quadrature().points();
            [qpPos, qpNeg] = rooftopMoM.basisFunctions.utils.generateQuadraturePointsUniformRooftop(e1, e2, shifts, pointsU, pointsV)

            % The first rooftop wedge is positive and has an origin in [0; 0; 0].
            rpPos = (1 / len) * (e * (e.' * (qpPos - shiftPos)))

            % The second one is negative and must be shifted by the total length of the rooftop (2 * len) in the e direction.
            rpNeg = (-1 / len) * (e * (e.' * (qpNeg - 2 * shiftNeg)))

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
