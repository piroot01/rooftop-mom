classdef GaussLegendre2RectQuadratureTest < unit.Base
% GaussLegendre2RectQuadratureTest: Tests of rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature.

    methods (Test)
        function testIntegratesPolynomialsExactly(testCase)
            order = 3;
            quadrature = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);
            [pointsU, pointsV] = quadrature.points();
            weights = quadrature.weights2D();

            for degreeU = 0:(2 * order - 1)
                for degreeV = 0:(2 * order - 1)
                    integral = (pointsU.^degreeU .* pointsV.^degreeV) * weights;
                    testCase.verifyEqual(integral, 1 / ((degreeU + 1) * (degreeV + 1)), 'RelTol', 1e-14);
                end
            end
        end

        function testUCoordinateRunsFastest(testCase)
            quadrature = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(2);
            [pointsU, pointsV] = quadrature.points();

            testCase.verifyEqual(pointsU(1:2), pointsU(3:4));
            testCase.verifyEqual(pointsV(1), pointsV(2));
            testCase.verifyLessThan(pointsU(1), pointsU(2));
        end
    end
end
