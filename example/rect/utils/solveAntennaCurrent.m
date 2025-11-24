function output = solveAntennaCurrent(solverContext, varargin)
    p = inputParser();

    addParameter(p, 'tolerance', 1e-12);
    addParameter(p, 'maxIterCount', 1000);
    addParameter(p, 'useSplitFFT', false);
    addParameter(p, 'genome', []);

    parse(p, varargin{:});

    tolerance = p.Results.tolerance;
    maxIterCount = p.Results.maxIterCount;
    useSplitFFT = p.Results.useSplitFFT;

    % Create the Green's matrix applicator
    if ~useSplitFFT
        greensMatrixApplicator = rooftopMoM.core.greensMatrix.FFTAcceleratedApplicator(solverContext.context, solverContext.greensMatrix);
    else
        greensMatrixApplicator = rooftopMoM.core.greensMatrix.SplitFFTAcceleratedApplicator(solverContext.context, solverContext.greensMatrix);
    end

    greensMatrixApplicator.generate();

    TRACE("solveAntennaCurrent: Green's matrix applicator generated.");

    % Instantiate the preconditioned conjugate gradient solver
    solver = rooftopMoM.core.solver.SpectrallyPreconditionedConjugateGradient(greensMatrixApplicator);

    % Finally solve for the current
    solverOutput = solver.solve( ...
        'excitation', solverContext.excitation, ...
        'tolerance', tolerance, ...
        'maxIterCount', maxIterCount ...
    );

    output = struct( ...
        'solverOutput', solverOutput ...
    );
end
