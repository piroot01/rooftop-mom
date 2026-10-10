classdef Resources
% Resources: User control of the memory and the threads of rooftopMoM. Before an object of the API allocates its data,
%   it checks that the data fits below the memory ceiling and into the free memory of the machine.
%       rooftopMoM.Resources.configure('memory', 2 * 2^30, 'threads', 8);
%       report = rooftopMoM.Resources.check(mesh, 'method', 'gmres', 'restart', 30);
%   rooftopMoM.Resources.reset() restores no ceiling and the automatic number of threads.

    methods (Static)
        function configure(options)
        % configure: Sets the memory ceiling and the number of computational threads.
        %
        %   Optional name-value inputs:
        %       'memory': Bytes the objects of rooftopMoM may hold together; the free memory of the machine still
        %           applies (default Inf, no ceiling).
        %       'threads': Number of computational threads of MATLAB, or 'automatic' (default 'automatic').

            arguments
                options.memory (1, 1) {mustBeReal, mustBeNonNan, mustBePositive} = Inf
                options.threads {mustBeThreadCount_} = 'automatic'
            end

            rooftopMoM.resources.Budget.instance().configure(double(options.memory));

            % FFTW keeps the thread count of its first call in the session as its maximum, so it is started with all
            % threads before the count is lowered.
            maxNumCompThreads('automatic');
            fftw('planner');
            maxNumCompThreads(options.threads);
        end

        function reset()
        % reset: Restores no memory ceiling and the automatic number of threads.

            rooftopMoM.Resources.configure();
        end

        function bytes = available()
        % available: Returns the free memory of the machine in bytes, NaN if it cannot be read.

            bytes = rooftopMoM.resources.readAvailableMemory();
        end

        function bytes = ceiling()
        % ceiling: Returns the memory ceiling in bytes, Inf for none.

            bytes = rooftopMoM.resources.Budget.instance().ceiling();
        end

        function bytes = held()
        % held: Returns the bytes held by the live objects of rooftopMoM.

            bytes = rooftopMoM.resources.Budget.instance().held();
        end

        function report = check(mesh, options)
        % check: Returns the memory a planned run needs, estimated before anything is built, and raises an error if
        %   its peak does not fit below the ceiling or into the free memory of the machine.
        %
        %   Input:
        %       mesh: rooftopMoM.Mesh.
        %
        %   Optional name-value inputs, as for Problem, Operator, Preconditioner and Solver:
        %       'chunkSize': Bytes of one chunk of the assembly (default as for Problem).
        %       'product': 'fft', 'splitFFT' or 'leanSplitFFT' (default 'splitFFT').
        %       'preconditioner': Logical flag, the run uses a circulant preconditioner (default true).
        %       'method': 'cg' or 'gmres' (default 'cg').
        %       'restart': Iterations between the restarts of 'gmres' (default 50).
        %       'maxIterationCount': Maximum number of iterations (default 1000).
        %
        %   Output:
        %       report: Struct of the bytes of the generators, the assembly chunk, the product and its working memory,
        %           the preconditioner and its working memory, the solver vectors and the peak of the run.

            arguments
                mesh (1, 1) rooftopMoM.Mesh
                options.chunkSize {mustBeScalarOrEmpty, mustBeReal, mustBePositive} = []
                options.product (1, :) char {mustBeMember(options.product, {'fft', 'splitFFT', 'leanSplitFFT'})} = ...
                    'splitFFT'
                options.preconditioner (1, 1) logical = true
                options.method (1, :) char {mustBeMember(options.method, {'cg', 'gmres'})} = 'cg'
                options.restart (1, 1) {mustBePositive} = 50
                options.maxIterationCount (1, 1) {mustBeInteger, mustBePositive} = 1000
            end

            % The memory does not depend on the quadrature.
            rooftops = rooftopMoM.basisFunctions.UniformRooftops(mesh.primitive(), ...
                rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(1));
            budget = rooftopMoM.resources.Budget.instance();

            report.generators = rooftopMoM.resources.estimateGenerators(rooftops);
            report.chunk = options.chunkSize;

            if isempty(report.chunk)
                report.chunk = budget.chooseChunk(report.generators);
            end

            [report.product, report.productWorking] = rooftopMoM.resources.estimateProduct(rooftops, options.product);
            [report.preconditioner, report.preconditionerWorking] = deal(0);

            if options.preconditioner
                [report.preconditioner, report.preconditionerWorking] = ...
                    rooftopMoM.resources.estimatePreconditioner(rooftops);
            end

            report.solver = rooftopMoM.resources.estimateSolver(rooftops.count(), options.method, ...
                options.preconditioner, options.restart, options.maxIterationCount);

            % The assembly and the solve do not overlap; within the solve one product runs at a time.
            report.peak = max(report.generators + report.chunk, report.generators + report.product ...
                + report.preconditioner + max(report.productWorking, report.preconditionerWorking) + report.solver);

            toText = @rooftopMoM.resources.utils.formatBytes;
            budget.require(report.peak, sprintf(['The run of %d unknowns (generators %s, chunk %s, product %s + %s ' ...
                'working, preconditioner %s + %s working, solver %s)'], rooftops.count(), toText(report.generators), ...
                toText(report.chunk), toText(report.product), toText(report.productWorking), ...
                toText(report.preconditioner), toText(report.preconditionerWorking), toText(report.solver)));
        end
    end
end

function mustBeThreadCount_(value)
    isAutomatic = (ischar(value) || isstring(value)) && strcmpi(value, 'automatic');
    isCount = isnumeric(value) && isscalar(value) && value >= 1 && value == round(value);

    if ~isAutomatic && ~isCount
        error('rooftopMoM:resources:invalidThreads', 'The threads are a positive integer or ''automatic''.');
    end
end
