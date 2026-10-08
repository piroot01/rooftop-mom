classdef ComputeRadarCrossSectionTest < unit.Base
% ComputeRadarCrossSectionTest: Tests of rooftopMoM.operator.computeRadarCrossSection.

    methods (Test)
        function testValue(testCase)
            % excitation.' * current = 2i.
            crossSection = rooftopMoM.operator.computeRadarCrossSection(2, [1; 1i], [1i; 1], 1);

            testCase.verifyEqual(crossSection, (2 * rooftopMoM.assembler.Constants.Z0)^2 * 4 / (4 * pi), ...
                'RelTol', 1e-15);
        end

        function testDoesNotDependOnTheIntensity(testCase)
            % The excitation and the current scale with the intensity.
            rng(3);
            excitation = complex(randn(5, 1), randn(5, 1));
            current = complex(randn(5, 1), randn(5, 1));

            weak = rooftopMoM.operator.computeRadarCrossSection(2 * pi, excitation, current, 1);
            strong = rooftopMoM.operator.computeRadarCrossSection(2 * pi, 3 * excitation, 3 * current, 3);

            testCase.verifyEqual(strong, weak, 'RelTol', 1e-14);
        end
    end
end
