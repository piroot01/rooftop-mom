classdef GenerateDeltaGapTest < unit.Base
% GenerateDeltaGapTest: Tests of rooftopMoM.operator.generateDeltaGap.

    methods (Test)
        function testVoltageTimesTheEdgeLength(testCase)
            rooftops = testCase.createRooftops(4, 3, 0.1, 0.2, 3);
            excitation = rooftopMoM.operator.generateDeltaGap(rooftops, [2, 5, 11], 3);

            % Rooftops 2 and 5 are x-directed with edges of 0.2, rooftop 11 is y-directed with an edge of 0.1.
            expected = zeros(17, 1);
            expected([2, 5]) = 3 * 0.2;
            expected(11) = 3 * 0.1;
            testCase.verifyEqual(excitation, expected);
        end
    end
end
