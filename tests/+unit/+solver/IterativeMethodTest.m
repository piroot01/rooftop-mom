classdef IterativeMethodTest < unit.solver.Base
% IterativeMethodTest: Tests of the behaviour every method shares through rooftopMoM.solver.IterativeMethod.

    properties (TestParameter)
        % Every method from a tolerance and an iteration limit. The preconditioner is the identity, and GMRES runs
        % without restarts, so every limit counts single iterations.
        createMethod = struct( ...
            'ConjugateGradient', @(tolerance, maxIterationCount) ...
                rooftopMoM.solver.ConjugateGradient(tolerance, maxIterationCount), ...
            'PreconditionedConjugateGradient', @(tolerance, maxIterationCount) ...
                rooftopMoM.solver.PreconditionedConjugateGradient(@(input) input, tolerance, maxIterationCount), ...
            'GeneralizedMinimumResidual', @(tolerance, maxIterationCount) ...
                rooftopMoM.solver.GeneralizedMinimumResidual(tolerance, maxIterationCount, Inf, []));
    end

    methods (Test)
        function testSolvesTheSystem(testCase, createMethod)
            [matrix, rightHandSide] = testCase.createSymmetricSystem(30);
            output = createMethod(1e-10, 100).solve(matrix, rightHandSide, []);
            expected = matrix \ rightHandSide;

            testCase.verifyTrue(output.converged());
            testCase.verifyLessThan(norm(output.solution() - expected) / norm(expected), 1e-9);
            testCase.verifyNumElements(output.residualNorms(), output.iterationCount() + 1);
        end

        function testOutputHoldsTheTrueResidual(testCase, createMethod)
            [matrix, rightHandSide] = testCase.createSymmetricSystem(30);
            output = createMethod(1e-6, 100).solve(matrix, rightHandSide, []);
            residualNorm = norm(rightHandSide - matrix * output.solution()) / norm(rightHandSide);

            testCase.verifyEqual(output.residualNorm(), residualNorm, 'RelTol', 1e-12);
            testCase.verifyLessThan(residualNorm, 1e-6);
        end

        function testZeroRightHandSideGivesTheZeroSolution(testCase, createMethod)
            output = createMethod(1e-10, 100).solve(eye(3), zeros(3, 1), ones(3, 1));

            testCase.verifyEqual(output.solution(), zeros(3, 1));
            testCase.verifyTrue(output.converged());
            testCase.verifyEqual([output.iterationCount(), output.residualNorm()], [0, 0]);
        end

        function testExactInitialGuessNeedsNoIteration(testCase, createMethod)
            [matrix, rightHandSide] = testCase.createSymmetricSystem(30);
            output = createMethod(1e-10, 100).solve(matrix, rightHandSide, matrix \ rightHandSide);

            testCase.verifyTrue(output.converged());
            testCase.verifyEqual(output.iterationCount(), 0);
        end

        function testAcceptsEveryOperatorForm(testCase, createMethod)
            [matrix, rightHandSide] = testCase.createSymmetricSystem(30);
            method = createMethod(1e-10, 100);
            expected = method.solve(matrix, rightHandSide, []).solution();

            testCase.verifyEqual(method.solve(@(input) matrix * input, rightHandSide, []).solution(), expected);
            testCase.verifyEqual(method.solve(struct('apply', @(input) matrix * input), rightHandSide, []).solution(), ...
                expected);
        end

        function testReportsNonConvergence(testCase, createMethod)
            sink = rooftopMoM.logging.sinks.MemorySink();
            rooftopMoM.logging.Logger.instance().configure(rooftopMoM.logging.LogLevel.WARNING, {sink});
            [matrix, rightHandSide] = testCase.createSymmetricSystem(30);
            method = createMethod(1e-14, 2);

            output = method.solve(matrix, rightHandSide, []);

            lines = sink.lines();
            testCase.verifyFalse(output.converged());
            testCase.verifyEqual(output.iterationCount(), 2);
            testCase.verifyNumElements(lines, 1);
            testCase.verifySubstring(lines{1}, ['[warning] [IterativeMethod.solve]: ', nameOf_(method), ...
                ' did not converge in 2 iterations']);
        end

        function testLogsTheStartAndTheResult(testCase, createMethod)
            sink = rooftopMoM.logging.sinks.MemorySink();
            rooftopMoM.logging.Logger.instance().configure(rooftopMoM.logging.LogLevel.DEBUG, {sink});
            [matrix, rightHandSide] = testCase.createSymmetricSystem(30);
            method = createMethod(1e-10, 100);

            method.solve(matrix, rightHandSide, []);

            lines = sink.lines();
            testCase.verifyNumElements(lines, 2);
            testCase.verifySubstring(lines{1}, ['[debug] [IterativeMethod.solve]: Solving 30 unknowns with ', ...
                nameOf_(method), ' to the relative residual 1.0e-10 in at most 100 iterations.']);
            testCase.verifySubstring(lines{2}, ['[info] [IterativeMethod.solve]: ', nameOf_(method), ' converged in ']);
        end
    end
end

function name = nameOf_(method)
    name = extractAfter(class(method), 'rooftopMoM.solver.');
end
