function output = solveAntennaCurrent(solverContext, varargin)
% solveAntennaCurrent: Solves the MoM system for the antenna current using the preconditioned conjugate gradient
%   method with a circulant preconditioner.
%
%   Inputs:
%       solverContext: Struct from prepareSolverContext with the fields context, greensMatrix and excitation.
%       varargin: Name-value pairs:
%           'tolerance': Relative residual tolerance (default 1e-12).
%           'maxIterCount': Maximum number of iterations (default 1000).
%           'useSplitFFT': Logical flag selecting the split-FFT greens matrix applicator (default false).
%           'sizeX': Face size in the x direction (default 1).
%           'sizeY': Face size in the y direction (default 1).
%           'preconditionerType': Circulant approximation, 'TChan' or 'Strang' (default 'Strang').
%
%   Output:
%       output: Struct with the field solverOutput holding the solver output.

    p = inputParser();

    addParameter(p, 'tolerance', 1e-12);
    addParameter(p, 'maxIterCount', 1000);
    addParameter(p, 'useSplitFFT', false);
    addParameter(p, 'sizeX', 1);
    addParameter(p, 'sizeY', 1);
    addParameter(p, 'preconditionerType', 'Strang');

    parse(p, varargin{:});

    tolerance = p.Results.tolerance;
    maxIterCount = p.Results.maxIterCount;

    % Create the Green's matrix applicator
    greensMatrixApplicator = generateGreensMatrixApplicator(solverContext, ...
        'sizeX', p.Results.sizeX, ...
        'sizeY', p.Results.sizeY, ...
        'useSplitFFT', p.Results.useSplitFFT ...
    );

    TRACE("solveAntennaCurrent: Green's matrix applicator generated.");

    % Create the spectral preconditioner applicator
    preconditionerApplicator = generatePreconditionerApplicator(solverContext, ...
        'sizeX', p.Results.sizeX, ...
        'sizeY', p.Results.sizeY, ...
        'preconditionerType', p.Results.preconditionerType ...
    );

    TRACE("solveAntennaCurrent: Preconditioner applicator generated.");

    % Instantiate the preconditioned conjugate gradient solver
    solver = rooftopMoM.core.solver.PreconditionedConjugateGradient(greensMatrixApplicator, preconditionerApplicator);

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
