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
            testCase.verifyEqual(problem.context().quadratureOrders(), [3, 3]);
        end

        function testQuadratureOrdersOfTheNearAndTheFarPairs(testCase)
            mesh = rooftopMoM.Mesh(4, 3, 0.1, 0.1);
            default = rooftopMoM.Problem(mesh, 2 * pi);
            chosen = rooftopMoM.Problem(mesh, 2 * pi, 'quadratureOrder', [12, 6]);

            % The excitation integrates with the far order.
            testCase.verifyEqual(default.context().quadratureOrders(), [16, 8]);
            testCase.verifyEqual(default.rooftops().quadrature().order(), 8);
            testCase.verifyEqual(chosen.context().quadratureOrders(), [12, 6]);
            testCase.verifyEqual(chosen.rooftops().quadrature().order(), 6);
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

            testCase.verifyTrue(any(contains(sink.lines(), ['[Problem.Problem]: 4 x 3 faces of 0.1 x 0.2 m, ' ...
                '17 rooftops, wavenumber 6.28319 rad/m, quadrature orders 3 near and 3 far.'])));
        end

        function testInvalidInputsRaiseErrors(testCase)
            mesh = rooftopMoM.Mesh(4, 3, 0.1, 0.1);

            testCase.verifyError(@() rooftopMoM.Problem(mesh, 0), 'MATLAB:validators:mustBePositive');
            testCase.verifyError(@() rooftopMoM.Problem(mesh, 2 * pi, 'quadratureOrder', 2.5), ...
                'MATLAB:validators:mustBeInteger');
            testCase.verifyError(@() rooftopMoM.Problem(mesh, 2 * pi, 'quadratureOrder', [16, 8, 4]), ...
                'rooftopMoM:problem:invalidQuadratureOrder');
            testCase.verifyError(@() rooftopMoM.Problem(mesh, 2 * pi, 'order', 3), ...
                'MATLAB:TooManyInputs');
        end
    end
end
