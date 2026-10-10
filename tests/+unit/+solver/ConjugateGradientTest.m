classdef ConjugateGradientTest < unit.solver.Base
% ConjugateGradientTest: Tests of rooftopMoM.solver.ConjugateGradient and its preconditioned variant.

    methods (Test)
        function testSolvesAComplexSymmetricSystem(testCase)
            [matrix, rightHandSide] = testCase.createSymmetricSystem(40);
            output = rooftopMoM.solver.ConjugateGradient(1e-10, 100).solve(matrix, rightHandSide, []);
            expected = matrix \ rightHandSide;

            testCase.verifyTrue(output.converged());
            testCase.verifyLessThan(norm(output.solution() - expected) / norm(expected), 1e-9);
        end

        function testStopsOnTheResidualNorm(testCase)
            [matrix, rightHandSide] = testCase.createSymmetricSystem(40);
            output = rooftopMoM.solver.ConjugateGradient(1e-8, 100).solve(matrix, rightHandSide, []);
            residualNorms = output.residualNorms();
            residualNorm = norm(rightHandSide - matrix * output.solution()) / norm(rightHandSide);

            testCase.verifyEqual(residualNorms(1), 1);
            testCase.verifyEqual(residualNorms(end), residualNorm, 'RelTol', 1e-4);
            testCase.verifyLessThan(residualNorm, 1e-8);
        end

        function testRecordsEveryIteration(testCase)
            [matrix, rightHandSide] = testCase.createSymmetricSystem(40);
            output = rooftopMoM.solver.ConjugateGradient(1e-8, 100).solve(matrix, rightHandSide, []);
            residualNorms = output.residualNorms();
            shorter = rooftopMoM.solver.ConjugateGradient(1e-8, output.iterationCount() - 1).solve(matrix, ...
                rightHandSide, []);

            testCase.verifyLessThan(residualNorms(end), 1e-8);
            testCase.verifyGreaterThanOrEqual(residualNorms(1:(end - 1)), 1e-8);
            testCase.verifyFalse(shorter.converged());
        end

        function testTracesEveryIteration(testCase)
            sink = rooftopMoM.logging.sinks.MemorySink();
            rooftopMoM.logging.Logger.instance().configure(rooftopMoM.logging.LogLevel.TRACE, {sink});
            [matrix, rightHandSide] = testCase.createSymmetricSystem(40);

            output = rooftopMoM.solver.ConjugateGradient(1e-8, 100).solve(matrix, rightHandSide, []);

            lines = sink.lines();
            testCase.verifyEqual(nnz(contains(lines, '[trace] [ConjugateGradient.iterate_]: Iteration ')), ...
                output.iterationCount());
        end

        function testPreconditionerCutsTheIterations(testCase)
            [matrix, rightHandSide] = testCase.createSymmetricSystem(100);

            % Scaling the rows and columns by 1 to 100 keeps the symmetry, but the diagonal preconditioner undoes it.
            scaling = logspace(0, 2, 100).';
            matrix = scaling .* matrix .* scaling.';
            diagonal = diag(matrix);

            plain = rooftopMoM.solver.ConjugateGradient(1e-8, 1000).solve(matrix, rightHandSide, []);
            preconditioned = rooftopMoM.solver.PreconditionedConjugateGradient(@(residual) residual ./ diagonal, ...
                1e-8, 1000).solve(matrix, rightHandSide, []);

            testCase.verifyTrue(preconditioned.converged());
            testCase.verifyLessThan(preconditioned.residualNorm(), 1e-8);
            testCase.verifyLessThan(preconditioned.iterationCount(), plain.iterationCount() / 2);
        end

        function testStopsOnABreakdown(testCase)
            % The first direction b has b.' * A * b = 0, so the step length is infinite.
            output = rooftopMoM.solver.ConjugateGradient(1e-8, 100).solve(diag([1, -1]), [1; 1], []);

            testCase.verifyFalse(output.converged());
            testCase.verifyEqual(output.iterationCount(), 1);
        end
    end
end
