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

        function [rooftopPoints, rooftopValues] = generateRooftopsDir(obj, dir)
        % generateRooftopsDirection: Generate rooftop basis for given axis direction.

            % Generate a single rooftop's geometry and values.
            [rooftopPoints, rooftopValues] = obj.generateRooftop_(dir);
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
