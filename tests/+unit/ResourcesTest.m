classdef ResourcesTest < unit.Base
% ResourcesTest: Tests of rooftopMoM.Resources and of the memory checks of the API.

    methods (Test)
        function testConfigureSetsTheCeilingAndTheThreads(testCase)
            rooftopMoM.Resources.configure('memory', 2^30, 'threads', 2);

            testCase.verifyEqual(rooftopMoM.Resources.ceiling(), 2^30);
            testCase.verifyEqual(maxNumCompThreads, 2);

            rooftopMoM.Resources.reset();

            testCase.verifyEqual(rooftopMoM.Resources.ceiling(), Inf);
        end

        function testProblemAboveTheCeilingRaisesError(testCase)
            rooftopMoM.Resources.configure('memory', 2^16);

            % The message names what the run needs, so the user knows what to give it.
            testCase.verifyError(@() rooftopMoM.Problem(rooftopMoM.Mesh(64, 64, 0.05, 0.05), 2 * pi), ...
                'rooftopMoM:resources:insufficientMemory');

            try
                rooftopMoM.Problem(rooftopMoM.Mesh(64, 64, 0.05, 0.05), 2 * pi);
            catch exception
                testCase.verifySubstring(exception.message, 'The assembly of 8064 rooftops (generators');
            end
        end

        function testOperatorAboveTheCeilingRaisesError(testCase)
            % The ceiling holds the assembly of the generators, but not the spectra of the product.
            mesh = rooftopMoM.Mesh(32, 32, 0.05, 0.05);
            rooftopMoM.Resources.configure('memory', testCase.estimateGenerators_(mesh) + 2^16);
            problem = rooftopMoM.Problem(mesh, 2 * pi, 'quadratureOrder', 2);

            testCase.verifyError(@() rooftopMoM.Operator(problem), 'rooftopMoM:resources:insufficientMemory');
        end

        function testDefaultChunkStaysBelowTheCeiling(testCase)
            mesh = rooftopMoM.Mesh(32, 32, 0.05, 0.05);
            rooftopMoM.Resources.configure('memory', testCase.estimateGenerators_(mesh) + 2^18);
            problem = rooftopMoM.Problem(mesh, 2 * pi, 'quadratureOrder', 2);

            testCase.verifyEqual(problem.context().chunkSize(), 2^18);
        end

        function testHeldMemoryFollowsTheLiveObjects(testCase)
            held = rooftopMoM.Resources.held();
            problem = rooftopMoM.Problem(rooftopMoM.Mesh(16, 16, 0.05, 0.05), 2 * pi, 'quadratureOrder', 2);
            operator = rooftopMoM.Operator(problem);

            testCase.verifyEqual(rooftopMoM.Resources.held(), held + problem.memory() + operator.memory());

            clear operator problem;

            testCase.verifyEqual(rooftopMoM.Resources.held(), held);
        end

        function testCheckMatchesTheBuiltRun(testCase)
            mesh = rooftopMoM.Mesh(24, 20, 0.05, 0.05);
            report = rooftopMoM.Resources.check(mesh, 'product', 'splitFFT');
            problem = rooftopMoM.Problem(mesh, 2 * pi, 'quadratureOrder', 2);
            operator = rooftopMoM.Operator(problem, 'product', 'splitFFT');
            preconditioner = rooftopMoM.Preconditioner(problem);

            testCase.verifyEqual([report.generators, report.chunk], [problem.memory(), problem.context().chunkSize()]);
            testCase.verifyEqual([report.product, report.productWorking], ...
                [operator.memory(), operator.workingMemory()]);
            testCase.verifyEqual([report.preconditioner, report.preconditionerWorking], ...
                [preconditioner.memory(), preconditioner.workingMemory()]);
        end

        function testCheckRaisesErrorAboveTheCeiling(testCase)
            % A chunk set by the user stays as it is, so a ceiling below the peak cannot be met.
            mesh = rooftopMoM.Mesh(64, 64, 0.05, 0.05);
            check = @() rooftopMoM.Resources.check(mesh, 'chunkSize', 2^20, 'method', 'gmres', 'restart', 20);
            report = check();
            rooftopMoM.Resources.configure('memory', report.peak - 1);

            testCase.verifyError(check, 'rooftopMoM:resources:insufficientMemory');
        end

        function testDefaultChunkAdaptsToTheCeiling(testCase)
            % The assembly with the default chunk is the peak of a small run; a lower ceiling only shrinks the chunk.
            mesh = rooftopMoM.Mesh(64, 64, 0.05, 0.05);
            report = rooftopMoM.Resources.check(mesh);
            rooftopMoM.Resources.configure('memory', report.peak - 2^20);

            adapted = rooftopMoM.Resources.check(mesh);

            testCase.verifyEqual(adapted.chunk, report.chunk - 2^20);
            testCase.verifyLessThanOrEqual(adapted.peak, report.peak - 2^20);
        end

        function testConfiguredChunkSizeHoldsForEveryProblem(testCase)
            mesh = rooftopMoM.Mesh(16, 16, 0.05, 0.05);
            rooftopMoM.Resources.configure('chunkSize', 2^20);

            problem = rooftopMoM.Problem(mesh, 2 * pi, 'quadratureOrder', 2);
            report = rooftopMoM.Resources.check(mesh);

            testCase.verifyEqual([problem.context().chunkSize(), report.chunk], [2^20, 2^20]);
        end

        function testChunkSizeOfAProblemOverridesTheConfiguredOne(testCase)
            rooftopMoM.Resources.configure('chunkSize', 2^20);

            problem = rooftopMoM.Problem(rooftopMoM.Mesh(16, 16, 0.05, 0.05), 2 * pi, 'quadratureOrder', 2, ...
                'chunkSize', 2^21);

            testCase.verifyEqual(problem.context().chunkSize(), 2^21);
        end

        function testConfiguredChunkSizeIsNotReducedToFitTheCeiling(testCase)
            % The chunk fits below the ceiling, but not next to the generators.
            mesh = rooftopMoM.Mesh(32, 32, 0.05, 0.05);
            rooftopMoM.Resources.configure('memory', testCase.estimateGenerators_(mesh) + 2^18, 'chunkSize', 2^18 + 1);

            testCase.verifyError(@() rooftopMoM.Problem(mesh, 2 * pi, 'quadratureOrder', 2), ...
                'rooftopMoM:resources:insufficientMemory');
        end

        function testResetRestoresTheAutomaticChunkSize(testCase)
            mesh = rooftopMoM.Mesh(16, 16, 0.05, 0.05);
            automatic = rooftopMoM.Resources.check(mesh).chunk;
            rooftopMoM.Resources.configure('chunkSize', 2^20);

            rooftopMoM.Resources.reset();

            testCase.verifyEqual(rooftopMoM.Resources.check(mesh).chunk, automatic);
        end

        function testInvalidInputsRaiseErrors(testCase)
            testCase.verifyError(@() rooftopMoM.Resources.configure('memory', -1), 'MATLAB:validators:mustBePositive');
            testCase.verifyError(@() rooftopMoM.Resources.configure('threads', 0), ...
                'rooftopMoM:resources:invalidThreads');
            testCase.verifyError(@() rooftopMoM.Resources.configure('threads', 'many'), ...
                'rooftopMoM:resources:invalidThreads');
            testCase.verifyError(@() rooftopMoM.Resources.configure('chunkSize', 0), ...
                'rooftopMoM:resources:invalidChunkSize');
            testCase.verifyError(@() rooftopMoM.Resources.configure('chunkSize', 'large'), ...
                'rooftopMoM:resources:invalidChunkSize');
            testCase.verifyError(@() rooftopMoM.Resources.configure('memory', 2^20, 'chunkSize', 2^21), ...
                'rooftopMoM:resources:invalidChunkSize');
        end
    end

    methods (Test, TestTags = {'Slow'})
        function testFFTUsesAllThreadsAgainAfterOne(testCase)
            % FFTW keeps the thread count of its first call in a MATLAB session as its maximum, so only a fresh session
            % shows it: one thread first, then all of them again.
            testCase.assumeGreaterThan(feature('numcores'), 2, 'Too few cores to tell one thread from all of them.');
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            command = sprintf(['addpath(''%s''); rooftopMoM.Resources.configure(''threads'', 1); fft2(ones(64)); ' ...
                'rooftopMoM.Resources.configure(); A = complex(randn(2048), randn(2048)); fft2(A); cpu = cputime; ' ...
                'wall = tic; for r = 1:20, fft2(A); end; fprintf(''busy cores %%.2f\\n'', (cputime - cpu) / toc(wall));'], ...
                fullfile(root, 'src'));

            [status, output] = system(sprintf('"%s" -batch "%s"', fullfile(matlabroot, 'bin', 'matlab'), command));
            busy = sscanf(regexp(output, 'busy cores [\d.]+', 'match', 'once'), 'busy cores %f');

            testCase.verifyEqual(status, 0, output);
            testCase.verifyGreaterThan(busy, 2, output);
        end
    end

    methods (Access = private)
        function bytes = estimateGenerators_(testCase, mesh)
            [countX, countY] = mesh.count();
            [sizeX, sizeY] = mesh.faceSize();
            bytes = rooftopMoM.resources.estimateGenerators(testCase.createRooftops(countX, countY, sizeX, sizeY, 2));
        end
    end
end
