classdef GenerateGaussLegendre1DTest < unit.Base
% GenerateGaussLegendre1DTest: Tests of rooftopMoM.basisFunctions.utils.generateGaussLegendre1D.

    properties (TestParameter)
        order = {1, 2, 5, 12};
    end

    methods (Test)
        function testIntegratesPolynomialsExactly(testCase, order)
            [points, weights] = rooftopMoM.basisFunctions.utils.generateGaussLegendre1D(order, 0, 1);

            for degree = 0:(2 * order - 1)
                testCase.verifyEqual(weights.' * points.^degree, 1 / (degree + 1), 'RelTol', 1e-14);
            end
        end

        function testPointsAreAscendingInsideTheInterval(testCase, order)
            points = rooftopMoM.basisFunctions.utils.generateGaussLegendre1D(order, -2, 3);

            testCase.verifySize(points, [order, 1]);
            testCase.verifyTrue(issorted(points));
            testCase.verifyTrue(all(points > -2 & points < 3));
        end

        function testKnownTwoPointRule(testCase)
            [points, weights] = rooftopMoM.basisFunctions.utils.generateGaussLegendre1D(2, -1, 1);

            testCase.verifyEqual(points, [-1; 1] / sqrt(3), 'AbsTol', 1e-15);
            testCase.verifyEqual(weights, [1; 1], 'AbsTol', 1e-15);
        end
    end
end
