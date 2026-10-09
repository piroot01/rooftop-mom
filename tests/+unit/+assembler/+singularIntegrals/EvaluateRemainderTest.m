classdef EvaluateRemainderTest < unit.Base
% EvaluateRemainderTest: Tests of rooftopMoM.assembler.singularIntegrals.evaluateRemainder.

    methods (Test)
        function testIsTheGreensFunctionWithoutTheExtractedPart(testCase)
            [distance, k] = deal([0.05, 0.3, 1, 2.5], 2 * pi);
            expected = exp(-1i * k * distance) ./ distance - (1 ./ distance - k^2 * distance / 2 + k^4 * distance.^3 / 24);

            testCase.verifyEqual(rooftopMoM.assembler.singularIntegrals.evaluateRemainder(distance, k), expected, ...
                'RelTol', 1e-13);
        end

        function testKeepsItsAccuracyAtSmallDistances(testCase)
            % The Taylor series -ik + i k^3 R^2 / 6 - i k^5 R^4 / 120 - k^6 R^5 / 720 + ..., which the direct
            % difference loses to cancellation.
            [distance, k] = deal([1e-8, 1e-5, 1e-3], 2 * pi);
            series = -1i * k + 1i * k^3 * distance.^2 / 6 - 1i * k^5 * distance.^4 / 120 - k^6 * distance.^5 / 720;

            testCase.verifyEqual(rooftopMoM.assembler.singularIntegrals.evaluateRemainder(distance, k), series, ...
                'RelTol', 1e-14);
        end
    end
end
