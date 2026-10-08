classdef ComputeInputImpedanceTest < unit.Base
% ComputeInputImpedanceTest: Tests of rooftopMoM.operator.computeInputImpedance.

    methods (Test)
        function testVoltageOverThePortCurrent(testCase)
            rooftops = testCase.createRooftops(4, 3, 0.1, 0.2, 3);
            current = zeros(17, 1);
            current([2, 5]) = [1 + 2i; 3];

            impedance = rooftopMoM.operator.computeInputImpedance(rooftops, current, [2, 5], 5);

            % The port edges are 0.2 long.
            testCase.verifyEqual(impedance, 5 / ((4 + 2i) * 0.2), 'RelTol', 1e-15);
        end

        function testDoesNotDependOnTheScale(testCase)
            % A strip dipole half a wavelength long, scaled together with the wavelength: the electrical problem and
            % so the input impedance stay the same.
            impedance = @(scale) testCase.solveDipole_(scale * 0.025, scale * 0.004, 2 * pi / scale);

            testCase.verifyEqual(impedance(10), impedance(1), 'RelTol', 1e-10);
        end
    end

    methods (Access = private)
        function impedance = solveDipole_(testCase, sizeX, sizeY, wavenumber)
            context = testCase.createContext(20, 1, sizeX, sizeY, 4, wavenumber);
            rooftops = context.rooftops();
            generators = rooftopMoM.assembler.greensMatrix.FullAssembler().generate(context);
            port = rooftopMoM.operator.utils.generateMiddlePort(rooftops, 'x');

            current = testCase.generateDenseMatrix(generators, rooftops) \ ...
                rooftopMoM.operator.generateDeltaGap(rooftops, port, 1);
            impedance = rooftopMoM.operator.computeInputImpedance(rooftops, current, port, 1);
        end
    end
end
