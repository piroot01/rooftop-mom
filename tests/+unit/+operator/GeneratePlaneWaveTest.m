classdef GeneratePlaneWaveTest < unit.Base
% GeneratePlaneWaveTest: Tests of rooftopMoM.operator.generatePlaneWave.

    methods (Test)
        function testNormalIncidenceGivesTheFaceArea(testCase)
            % A rooftop of unit height covers half of each of its two faces, so the area of one face in total.
            rooftops = testCase.createRooftops(4, 3, 0.1, 0.2, 3);
            excitation = rooftopMoM.operator.generatePlaneWave(rooftops, 2 * pi, [0; 0; -1], [0.6; 0.8; 0], 2);

            testCase.verifyEqual(excitation, [repmat(2 * 0.6 * 0.02, 9, 1); repmat(2 * 0.8 * 0.02, 8, 1)], ...
                'RelTol', 1e-14);
        end

        function testStripHasOnlyItsDirection(testCase)
            rooftops = testCase.createRooftops(6, 1, 0.1, 0.02, 3);
            excitation = rooftopMoM.operator.generatePlaneWave(rooftops, 2 * pi, [0; 0; -1], [1; 0; 0], 1);

            testCase.verifyEqual(excitation, repmat(0.1 * 0.02, 5, 1), 'RelTol', 1e-14);
        end

        function testPhaseFollowsTheDirectionOfPropagation(testCase)
            % Neighbouring rooftops differ by the phase of the wave over one face, exp(-1i * k * direction.' * shift).
            rooftops = testCase.createRooftops(4, 3, 0.1, 0.2, 3);
            direction = [1; 2; -2] / 3;
            excitation = rooftopMoM.operator.generatePlaneWave(rooftops, 2 * pi, direction, [2; -1; 0] / sqrt(5), 1);
            partX = reshape(excitation(1:9), 3, 3);
            partY = reshape(excitation(10:17), 4, 2);

            testCase.verifyEqual(partX(2:end, :) ./ partX(1:(end - 1), :), ...
                repmat(exp(-1i * 2 * pi * direction(1) * 0.1), 2, 3), 'RelTol', 1e-13);
            testCase.verifyEqual(partY(:, 2) ./ partY(:, 1), repmat(exp(-1i * 2 * pi * direction(2) * 0.2), 4, 1), ...
                'RelTol', 1e-13);
        end

        function testMatchesAdaptiveIntegration(testCase)
            % One x-directed rooftop, rising over the face [0, 0.1] x [0, 0.2] and falling over the next one.
            rooftops = testCase.createRooftops(2, 1, 0.1, 0.2, 8);
            direction = [1; 2; -2] / 3;
            polarization = [2; -1; 0] / sqrt(5);
            excitation = rooftopMoM.operator.generatePlaneWave(rooftops, 2 * pi, direction, polarization, 1);

            wave = @(x, y) polarization(1) * exp(-1i * 2 * pi * (direction(1) * x + direction(2) * y));
            rising = integral2(@(x, y) x / 0.1 .* wave(x, y), 0, 0.1, 0, 0.2, 'AbsTol', 1e-16, 'RelTol', 1e-13);
            falling = integral2(@(x, y) (0.2 - x) / 0.1 .* wave(x, y), 0.1, 0.2, 0, 0.2, 'AbsTol', 1e-16, ...
                'RelTol', 1e-13);

            testCase.verifyEqual(excitation, rising + falling, 'RelTol', 1e-12);
        end
    end
end
