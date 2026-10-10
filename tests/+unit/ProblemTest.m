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
                '17 rooftops, wavenumber 6.28319 rad/m, ka 2.265, 5 faces per wavelength, quadrature orders 3 near ' ...
                'and 3 far.'])));
        end

        function testCoarseMeshLogsAWarning(testCase)
            sink = rooftopMoM.logging.sinks.MemorySink();
            rooftopMoM.logging.Logger.instance().configure(rooftopMoM.logging.LogLevel.WARNING, {sink});

            % Faces of 0.2 m give five per wavelength of 1 m, faces of 0.1 m exactly ten.
            rooftopMoM.Problem(rooftopMoM.Mesh(4, 3, 0.1, 0.2), 2 * pi, 'quadratureOrder', 3);
            rooftopMoM.Problem(rooftopMoM.Mesh(4, 3, 0.1, 0.1), 2 * pi, 'quadratureOrder', 3);

            lines = sink.lines();
            testCase.verifyNumElements(lines, 1);
            testCase.verifySubstring(lines{1}, ['[warning] [Problem.Problem]: The mesh has 5 faces per wavelength, ' ...
                'fewer than 10; at this wavenumber the faces should be at most 0.1 m.']);
        end

        function testElectricalSizeSetsTheWavenumber(testCase)
            mesh = rooftopMoM.Mesh(4, 3, 0.1, 0.2);

            problem = rooftopMoM.Problem(mesh, 'electricalSize', 0.5, 'quadratureOrder', 3);

            testCase.verifyEqual(problem.wavenumber(), 0.5 / mesh.radius(), 'RelTol', 1e-15);
            testCase.verifyEqual(problem.electricalSize(), 0.5, 'RelTol', 1e-15);
        end

        function testDensityCountsTheLongerFaceSide(testCase)
            % Faces of 0.08 m along y, 12.5 per wavelength of 1 m.
            problem = rooftopMoM.Problem(rooftopMoM.Mesh(4, 3, 0.05, 0.08), 2 * pi, 'quadratureOrder', 3);

            testCase.verifyEqual(problem.density(), 12.5, 'RelTol', 1e-15);
        end

        function testFromElectricalSizeMeshesThePlate(testCase)
            % ka = 4 at the wavelength 1 m: the diagonal is 4 / pi, the sides 1.139 and 0.569 m, which take 13.7 and
            % 6.8 faces at 12 per wavelength.
            problem = rooftopMoM.Problem.fromElectricalSize(4, 'density', 12, 'aspect', 2, 'quadratureOrder', 3);
            [countX, countY] = problem.mesh().count();
            [sizeX, sizeY] = problem.mesh().faceSize();

            testCase.verifyEqual([countX, countY], [14, 7]);
            testCase.verifyEqual(countX * sizeX / (countY * sizeY), 2, 'RelTol', 1e-14);
            testCase.verifyEqual(problem.wavenumber(), 2 * pi);
            testCase.verifyEqual(problem.electricalSize(), 4, 'RelTol', 1e-14);
            testCase.verifyGreaterThanOrEqual(problem.density(), 12);
        end

        function testFromElectricalSizeRoundsTheCountsToFastFFTSizes(testCase)
            % ka = 4.6 needs 10.35 faces per side at 10 per wavelength; the 11 faces that reach it round up to 12.
            problem = rooftopMoM.Problem.fromElectricalSize(4.6, 'density', 10, 'quadratureOrder', 3);
            [countX, countY] = problem.mesh().count();

            testCase.verifyEqual([countX, countY], [12, 12]);
            testCase.verifyGreaterThanOrEqual(problem.density(), 10);
        end

        function testFromElectricalSizeMeshesAStripWithOneFaceAcross(testCase)
            % A half-wave strip 250 times as long as wide: ka = pi / 2, 10 faces along it and 1 across.
            problem = rooftopMoM.Problem.fromElectricalSize(pi / 2, 'density', 20, 'aspect', 250, 'quadratureOrder', 3);
            [countX, countY] = problem.mesh().count();

            testCase.verifyEqual([countX, countY], [10, 1]);
        end

        function testFromElectricalSizeRefusesASingleFace(testCase)
            testCase.verifyError(@() rooftopMoM.Problem.fromElectricalSize(0.1), 'rooftopMoM:problem:singleFace');
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
            testCase.verifyError(@() rooftopMoM.Problem(mesh), 'rooftopMoM:problem:invalidWavenumber');
            testCase.verifyError(@() rooftopMoM.Problem(mesh, 2 * pi, 'electricalSize', 2), ...
                'rooftopMoM:problem:invalidWavenumber');
            testCase.verifyError(@() rooftopMoM.Problem(mesh, 'electricalSize', -1), ...
                'MATLAB:validators:mustBePositive');
            testCase.verifyError(@() rooftopMoM.Problem.fromElectricalSize(4, 'density', 0), ...
                'MATLAB:validators:mustBePositive');
        end
    end
end
