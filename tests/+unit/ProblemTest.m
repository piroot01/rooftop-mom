classdef ProblemTest < unit.Base
% ProblemTest: Tests of rooftopMoM.Problem.

    methods (Test)
        function testGetters(testCase)
            mesh = rooftopMoM.Mesh(4, 3, 0.1, 0.1);
            problem = rooftopMoM.Problem(mesh, 2 * pi, 'quadratureOrder', 3);

            testCase.verifyEqual(problem.mesh(), mesh);
            testCase.verifyEqual(problem.wavenumber(), 2 * pi);
            testCase.verifyEqual(problem.rooftopCount(), 17);
            testCase.verifyEqual(problem.rooftops().quadrature().order(), 3);
        end

        function testGeneratorsAreAssembledOnce(testCase)
            sink = rooftopMoM.logging.sinks.MemorySink();
            rooftopMoM.logging.Logger.instance().configure(rooftopMoM.logging.LogLevel.INFO, {sink});
            problem = rooftopMoM.Problem(rooftopMoM.Mesh(4, 3, 0.1, 0.1), 2 * pi, 'quadratureOrder', 3);

            testCase.verifyEqual(problem.memory(), 0);

            generators = problem.generators();

            testCase.verifyEqual(problem.generators(), generators);
            testCase.verifyEqual(problem.memory(), generators.memory());
            testCase.verifyNumElements(sink.lines(), 1);
        end

        function testLogsTheProblem(testCase)
            sink = rooftopMoM.logging.sinks.MemorySink();
            rooftopMoM.logging.Logger.instance().configure(rooftopMoM.logging.LogLevel.DEBUG, {sink});

            rooftopMoM.Problem(rooftopMoM.Mesh(4, 3, 0.1, 0.2), 2 * pi, 'quadratureOrder', 3);

            lines = sink.lines();
            testCase.verifySubstring(lines{1}, ['[Problem.Problem]: 4 x 3 faces of 0.1 x 0.2 m, 17 rooftops, ' ...
                'wavenumber 6.28319 rad/m, quadrature order 3.']);
        end

        function testInvalidInputsRaiseErrors(testCase)
            mesh = rooftopMoM.Mesh(4, 3, 0.1, 0.1);

            testCase.verifyError(@() rooftopMoM.Problem(mesh, 0), 'MATLAB:validators:mustBePositive');
            testCase.verifyError(@() rooftopMoM.Problem(mesh, 2 * pi, 'quadratureOrder', 2.5), ...
                'MATLAB:validators:mustBeInteger');
            testCase.verifyError(@() rooftopMoM.Problem(mesh, 2 * pi, 'order', 3), ...
                'MATLAB:TooManyInputs');
        end
    end
end
