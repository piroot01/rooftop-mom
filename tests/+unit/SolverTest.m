classdef SolverTest < unit.solver.Base
% SolverTest: Tests of rooftopMoM.Solver.

    properties (TestParameter)
        method = {'cg', 'gmres'};
    end

    methods (Test)
        function testSolvesAPlaneWaveOnAPlate(testCase, method)
            problem = testCase.problem_();
            excitation = rooftopMoM.Excitation.planeWave(problem);
            result = rooftopMoM.Solver('method', method, 'tolerance', 1e-10).solve(rooftopMoM.Operator(problem), ...
                excitation);
            matrix = testCase.generateDenseMatrix(problem.generators(), problem.rooftops());

            % The dense matrix checks the residual independently of the FFT product.
            testCase.verifyTrue(result.converged());
            testCase.verifyLessThan(norm(matrix * result.solution() - excitation.vector()) / ...
                norm(excitation.vector()), 1e-10);
            testCase.verifyEqual(result.excitation(), excitation);
        end

        function testPreconditionerCutsTheIterations(testCase, method)
            problem = testCase.problem_();
            operator = rooftopMoM.Operator(problem);
            excitation = rooftopMoM.Excitation.planeWave(problem);
            solver = rooftopMoM.Solver('method', method, 'tolerance', 1e-8);

            plain = solver.solve(operator, excitation);
            preconditioned = solver.solve(operator, excitation, 'preconditioner', rooftopMoM.Preconditioner(problem));

            testCase.verifyTrue(preconditioned.converged());
            testCase.verifyLessThan(preconditioned.iterationCount(), plain.iterationCount());
        end

        function testSolvesAPlainSystem(testCase, method)
            % Any matrix works, without a Problem.
            [matrix, rightHandSide] = testCase.createSymmetricSystem(20);
            result = rooftopMoM.Solver('method', method, 'tolerance', 1e-10).solve(matrix, rightHandSide);

            testCase.verifyTrue(result.converged());
            testCase.verifyLessThan(norm(matrix * result.solution() - rightHandSide) / norm(rightHandSide), 1e-10);
            testCase.verifyEmpty(result.excitation());
        end

        function testChecksTheSizes(testCase)
            solver = rooftopMoM.Solver();

            testCase.verifyError(@() solver.solve(eye(3), ones(4, 1)), 'rooftopMoM:solver:sizeMismatch');
            testCase.verifyError(@() solver.solve(rooftopMoM.Operator(testCase.problem_()), ones(4, 1)), ...
                'rooftopMoM:solver:sizeMismatch');
            testCase.verifyError(@() solver.solve(eye(3), ones(3, 1), 'preconditioner', eye(4)), ...
                'rooftopMoM:solver:sizeMismatch');
            testCase.verifyError(@() solver.solve(eye(3), ones(3, 1), 'initialGuess', ones(4, 1)), ...
                'rooftopMoM:solver:sizeMismatch');
        end

        function testInvalidInputsRaiseErrors(testCase)
            solver = rooftopMoM.Solver();

            testCase.verifyError(@() solver.solve('matrix', ones(3, 1)), 'rooftopMoM:solver:invalidOperator');
            testCase.verifyError(@() solver.solve(eye(3), {1, 2, 3}), 'rooftopMoM:solver:invalidRightHandSide');
            testCase.verifyError(@() solver.solve(eye(3), [1; NaN; 3]), 'rooftopMoM:solver:invalidRightHandSide');
            testCase.verifyError(@() rooftopMoM.Solver('method', 'bicg'), 'MATLAB:validators:mustBeMember');
            testCase.verifyError(@() rooftopMoM.Solver('tolerance', 2), 'MATLAB:validators:mustBeLessThan');
            testCase.verifyError(@() rooftopMoM.Solver('restart', 2.5), 'rooftopMoM:solver:invalidRestart');
        end
    end

    methods (Access = private)
        function problem = problem_(~)
            problem = rooftopMoM.Problem(rooftopMoM.Mesh(8, 8, 0.05, 0.05), 2 * pi, 'quadratureOrder', 3);
        end
    end
end
