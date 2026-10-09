classdef GenerateExtractedMomentsTest < unit.Base
% GenerateExtractedMomentsTest: Tests of rooftopMoM.assembler.singularIntegrals.generateExtractedMoments.

    methods (Test)
        function testCombinesThePowerMoments(testCase)
            [a, b, k] = deal(0.1, 0.15, 2 * pi);
            moments = @(power) rooftopMoM.assembler.singularIntegrals.generatePowerMoments(a, b, power);

            testCase.verifyEqual(rooftopMoM.assembler.singularIntegrals.generateExtractedMoments(a, b, k), ...
                moments(-1) - k^2 / 2 * moments(1) + k^4 / 24 * moments(3), 'RelTol', 1e-15);
        end

        function testStaticLimitIsTheMomentOfOneOverR(testCase)
            testCase.verifyEqual(rooftopMoM.assembler.singularIntegrals.generateExtractedMoments(0.1, 0.15, 0), ...
                rooftopMoM.assembler.singularIntegrals.generatePowerMoments(0.1, 0.15, -1));
        end
    end
end
