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
        end

        function testRadarCrossSectionOfAPlaneWave(testCase)
            problem = rooftopMoM.Problem(rooftopMoM.Mesh(6, 6, 0.05, 0.05), 2 * pi, 'quadratureOrder', 3);
            excitation = rooftopMoM.Excitation.planeWave(problem, 'intensity', 3);
            result = rooftopMoM.Solver('tolerance', 1e-10).solve(rooftopMoM.Operator(problem), excitation);

            testCase.verifyEqual(result.radarCrossSection(), rooftopMoM.operator.computeRadarCrossSection(2 * pi, ...
                excitation.vector(), result.solution(), 3));
            testCase.verifyError(@() result.inputImpedance(), 'rooftopMoM:result:wrongExcitation');
        end

        function testPlainSystemHasNoSourceQuantities(testCase)
            result = rooftopMoM.Solver().solve(2 * eye(3), ones(3, 1));

            testCase.verifyEqual(result.solution(), 0.5 * ones(3, 1), 'RelTol', 1e-14);
            testCase.verifyError(@() result.inputImpedance(), 'rooftopMoM:result:wrongExcitation');
        end
    end
end
