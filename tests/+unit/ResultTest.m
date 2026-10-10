classdef ResultTest < unit.Base
% ResultTest: Tests of rooftopMoM.Result.

    methods (Test)
        function testInputImpedanceOfADeltaGap(testCase)
            problem = rooftopMoM.Problem(rooftopMoM.Mesh(20, 1, 0.025, 0.004), 2 * pi, 'quadratureOrder', 3);
            excitation = rooftopMoM.Excitation.deltaGap(problem);
            result = rooftopMoM.Solver('tolerance', 1e-10).solve(rooftopMoM.Operator(problem), excitation);
            parameters = excitation.parameters();

            testCase.verifyEqual(result.inputImpedance(), rooftopMoM.operator.computeInputImpedance( ...
                problem.rooftops(), result.solution(), parameters.port, parameters.voltage));
            testCase.verifyError(@() result.radarCrossSection(), 'rooftopMoM:result:wrongExcitation');
            testCase.verifyError(@() result.normalizedRadarCrossSection(), 'rooftopMoM:result:wrongExcitation');
        end

        function testRadarCrossSectionOfAPlaneWave(testCase)
            problem = rooftopMoM.Problem(rooftopMoM.Mesh(6, 6, 0.05, 0.05), 2 * pi, 'quadratureOrder', 3);
            excitation = rooftopMoM.Excitation.planeWave(problem, 'intensity', 3);
            result = rooftopMoM.Solver('tolerance', 1e-10).solve(rooftopMoM.Operator(problem), excitation);

            testCase.verifyEqual(result.radarCrossSection(), rooftopMoM.operator.computeRadarCrossSection(2 * pi, ...
                excitation.vector(), result.solution(), 3));
            testCase.verifyError(@() result.inputImpedance(), 'rooftopMoM:result:wrongExcitation');
        end

        function testNormalizedRadarCrossSectionDoesNotDependOnTheScale(testCase)
            % The same plate at the wavelengths 1 m and 2 m: the cross section grows with the squared wavelength, and
            % its ratio to the squared wavelength stays.
            small = solvePlate_(1);
            large = solvePlate_(2);

            testCase.verifyEqual(large.radarCrossSection(), 4 * small.radarCrossSection(), 'RelTol', 1e-8);
            testCase.verifyEqual(large.normalizedRadarCrossSection(), small.normalizedRadarCrossSection(), ...
                'RelTol', 1e-8);
            testCase.verifyEqual(small.normalizedRadarCrossSection(), small.radarCrossSection(), 'RelTol', 1e-15);
        end

        function testPlainSystemHasNoSourceQuantities(testCase)
            result = rooftopMoM.Solver().solve(2 * eye(3), ones(3, 1));

            testCase.verifyEqual(result.solution(), 0.5 * ones(3, 1), 'RelTol', 1e-14);
            testCase.verifyError(@() result.inputImpedance(), 'rooftopMoM:result:wrongExcitation');
        end
    end
end

function result = solvePlate_(scale)
% solvePlate_: Returns the result of a plane wave on a plate of 6 x 6 faces at the wavelength scale, each face a
%   tenth of the wavelength wide.

    problem = rooftopMoM.Problem(rooftopMoM.Mesh(6, 6, 0.1 * scale, 0.1 * scale), 2 * pi / scale, ...
        'quadratureOrder', 3);
    result = rooftopMoM.Solver('tolerance', 1e-10).solve(rooftopMoM.Operator(problem), ...
        rooftopMoM.Excitation.planeWave(problem));
end
