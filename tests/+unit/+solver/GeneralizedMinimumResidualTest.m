classdef GeneralizedMinimumResidualTest < unit.solver.Base
% GeneralizedMinimumResidualTest: Tests of rooftopMoM.solver.GeneralizedMinimumResidual.

    methods (Test)
        function testSolvesANonsymmetricSystem(testCase)
            [matrix, rightHandSide] = testCase.createNonsymmetricSystem(40);
            output = rooftopMoM.solver.GeneralizedMinimumResidual(1e-10, 100, 10, []).solve(matrix, rightHandSide, []);
            expected = matrix \ rightHandSide;

            testCase.verifyTrue(output.converged());
            testCase.verifyGreaterThan(output.iterationCount(), 10);
            testCase.verifyLessThan(norm(output.solution() - expected) / norm(expected), 1e-9);
        end

        function testStopsAtTheIterationLimit(testCase)
            [matrix, rightHandSide] = testCase.createNonsymmetricSystem(60);
            output = rooftopMoM.solver.GeneralizedMinimumResidual(1e-14, 7, 5, []).solve(matrix, rightHandSide, []);

            testCase.verifyFalse(output.converged());
            testCase.verifyEqual(output.iterationCount(), 7);
        end

        function testRestartBeyondTheSystemSizeIsCapped(testCase)
            % A basis never needs more vectors than unknowns.
            [matrix, rightHandSide] = testCase.createNonsymmetricSystem(20);
            solve = @(restart) rooftopMoM.solver.GeneralizedMinimumResidual(1e-10, 100, restart, []).solve(matrix, ...
                rightHandSide, []);

            equal = solve(20);
            longer = solve(50);

            testCase.verifyTrue(equal.converged());
            testCase.verifyEqual(equal.solution(), longer.solution());
        end

        function testResidualIsThatOfThePreconditionedSystem(testCase)
            [matrix, rightHandSide] = testCase.createNonsymmetricSystem(60);

            % Scaling the rows by 1 to 1000 spoils the condition, but the diagonal preconditioner undoes it.
            matrix = logspace(0, 3, 60).' .* matrix;
            precondition = @(residual) residual ./ diag(matrix);

            plain = rooftopMoM.solver.GeneralizedMinimumResidual(1e-10, 200, 20, []).solve(matrix, rightHandSide, []);
            output = rooftopMoM.solver.GeneralizedMinimumResidual(1e-10, 200, 20, precondition).solve(matrix, ...
                rightHandSide, []);
            residualNorms = output.residualNorms();
            expected = norm(precondition(rightHandSide - matrix * output.solution())) / norm(precondition(rightHandSide));

            testCase.verifyTrue(output.converged());
            testCase.verifyEqual(residualNorms(1), 1, 'AbsTol', 1e-14);
            testCase.verifyEqual(residualNorms(end), expected, 'AbsTol', 1e-12);
            testCase.verifyLessThan(output.iterationCount(), plain.iterationCount());
        end

        function testTracesEveryIteration(testCase)
            sink = rooftopMoM.logging.sinks.MemorySink();
            rooftopMoM.logging.Logger.instance().configure(rooftopMoM.logging.LogLevel.TRACE, {sink});
            [matrix, rightHandSide] = testCase.createNonsymmetricSystem(40);

            output = rooftopMoM.solver.GeneralizedMinimumResidual(1e-10, 100, 10, []).solve(matrix, rightHandSide, []);

            lines = sink.lines();
            testCase.verifyEqual(nnz(contains(lines, ...
                '[trace] [GeneralizedMinimumResidual.iterate_]: Iteration ')), output.iterationCount());
        end
    end
end
