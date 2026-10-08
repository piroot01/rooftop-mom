classdef ExcitationTest < unit.Base
% ExcitationTest: Tests of rooftopMoM.Excitation.

    methods (Test)
        function testPlaneWaveNormalizesItsVectors(testCase)
            problem = testCase.problem_(4, 3);
            excitation = rooftopMoM.Excitation.planeWave(problem, 'direction', [0, 0, -2], ...
                'polarization', [0; 3; 0], 'intensity', 2);
            parameters = excitation.parameters();

            testCase.verifyEqual(excitation.type(), 'planeWave');
            testCase.verifyEqual([parameters.direction, parameters.polarization], [0, 0; 0, 1; -1, 0]);
            testCase.verifyEqual(excitation.vector(), rooftopMoM.operator.generatePlaneWave(problem.rooftops(), ...
                2 * pi, [0; 0; -1], [0; 1; 0], 2));
        end

        function testPlaneWaveNeedsPerpendicularVectors(testCase)
            problem = testCase.problem_(4, 3);

            testCase.verifyError(@() rooftopMoM.Excitation.planeWave(problem, 'direction', [1; 0; 0], ...
                'polarization', [1; 1; 0]), 'rooftopMoM:excitation:notPerpendicular');
            testCase.verifyError(@() rooftopMoM.Excitation.planeWave(problem, 'polarization', [0; 0; 0]), ...
                'rooftopMoM:excitation:zeroVector');
        end

        function testPlaneWaveAcceptsARotatedFrame(testCase)
            % The axes of a frame built from three rotations are perpendicular up to the rounding of the products.
            angles = [0.3, 1.1, 2.5];
            rotationX = [1, 0, 0; 0, cos(angles(1)), -sin(angles(1)); 0, sin(angles(1)), cos(angles(1))];
            rotationY = [cos(angles(2)), 0, sin(angles(2)); 0, 1, 0; -sin(angles(2)), 0, cos(angles(2))];
            rotationZ = [cos(angles(3)), -sin(angles(3)), 0; sin(angles(3)), cos(angles(3)), 0; 0, 0, 1];
            frame = rotationZ * rotationY * rotationX;

            excitation = rooftopMoM.Excitation.planeWave(testCase.problem_(4, 3), 'direction', frame(:, 3), ...
                'polarization', frame(:, 1));

            testCase.verifyEqual(excitation.type(), 'planeWave');
        end

        function testDeltaGapDefaultsToTheMiddlePort(testCase)
            problem = testCase.problem_(6, 1);
            excitation = rooftopMoM.Excitation.deltaGap(problem, 'voltage', 2);
            parameters = excitation.parameters();

            testCase.verifyEqual(excitation.type(), 'deltaGap');
            testCase.verifyEqual([parameters.port, parameters.voltage], [3, 2]);
            testCase.verifyEqual(excitation.vector(), rooftopMoM.operator.generateDeltaGap(problem.rooftops(), 3, 2));
        end

        function testDeltaGapChecksThePort(testCase)
            problem = testCase.problem_(6, 1);

            testCase.verifyError(@() rooftopMoM.Excitation.deltaGap(problem, 'port', [2, 2]), ...
                'rooftopMoM:excitation:invalidPort');
            testCase.verifyError(@() rooftopMoM.Excitation.deltaGap(problem, 'port', 6), ...
                'rooftopMoM:excitation:invalidPort');
            testCase.verifyError(@() rooftopMoM.Excitation.deltaGap(problem, 'direction', 'y'), ...
                'rooftopMoM:excitation:noRooftops');
        end
    end

    methods (Access = private)
        function problem = problem_(~, countX, countY)
            problem = rooftopMoM.Problem(rooftopMoM.Mesh(countX, countY, 0.1, 0.1), 2 * pi, 'quadratureOrder', 3);
        end
    end
end
