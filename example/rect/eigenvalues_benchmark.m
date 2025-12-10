%% eigenvalues benchmark
%
%   Description:

clc; clearvars; close all;

% Add required paths.
currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', '..', 'src');
addpath(srcDir);

utilsDir = fullfile(currentDir, 'utils');
addpath(utilsDir);

loggerDir = fullfile(currentDir, '..', '..', '3rd', 'matlab-logger', 'src');
addpath(loggerDir);

% Configure the logger.
logger.Logger.configure( ...
    'timeFormat', 'ss.SSS', ...
    'minLevel', logger.LogLevel.TRACE, ...
    'sinks', { ...
        logger.sink.ConsoleSink(), ...
    } ...
);

%% Prepare benchmark inputs
electricalSize = 2;

% Number of requested eigenvalues
eigCount = 10;
whichEigs = 'smallestabs';

% Constants used to compute the wavenumber and antenna size.
% Number of rooftops per wavelength to aim for.
rooftopCountPerOneWavelength = 10;

% Antennas with electrical size below this threshold will use the same
% discretization as an antenna whose electrical size equals the threshold.
electricalSizeThreshold = 5;

% Antenna side aspect ratio (sizeX:sizeY) will be 1:antennaSidesRatio.
antennaSidesRatio = 2;

% Order of the Gauss–Legendre quadrature.
quadratureOrder = 3;

% Define the incident plane wave.
direction = [0; 0; -1];
polarization = [0; 1; 0];
waveIntensity = 1;

% Memory allocation for a single chunk (number of elements, complex doubles).
chunkSize = 5e6;

% Iterative solver settings.
tolerance = 1e-12;
maxIterCount = 1e4;

%% Benchmark body
[countX, countY, sizeX, sizeY] = generateAntennaDimensions(electricalSize, antennaSidesRatio, rooftopCountPerOneWavelength);

% Compute the wavenumber for the current antenna geometry.
wavenumber = 2 * electricalSize / sqrt(sizeX^2 * countX^2 + sizeY^2 * countY^2);

DEBUG('RCS_benchmark: Preparing shared data for the solvers.');

solverContext = prepareSolverContext( ...
    'countX', countX, ...
    'countY', countY, ...
    'sizeX', sizeX, ...
    'sizeY', sizeY, ...
    'quadratureOrder', quadratureOrder, ...
    'wavenumber', wavenumber, ...
    'chunkSize', chunkSize, ...
    'waveDirection', direction, ...
    'wavePolarization', polarization, ...
    'waveIntensity', waveIntensity ...
);

DEBUG('RCS_benchmark: Computing with the FFT preconditioned CG solver');

applicators = generateApplicators(solverContext, ...
    'sizeX', sizeX, ...
    'sizeY', sizeY, ...
    'preconditionerType', 'TChan' ...
);

%% Test solver


%% First eigenvalue problem (eig1): imag(Z) * I = lambda * real(Z) * I
[lambda, V, info] = solve_eig1(eigCount, ...
    applicators, ...
    struct( ...
        'innerTolerance', tolerance, ...
        'innerMaxIterCount', maxIterCount, ...
        'whichEigs', whichEigs ...
    ) ...
);

disp(lambda);

function [lambda, V, info] = solve_eig1(eigCount, applicators, options)
% solve_eig1: Solves the imag(Z) * x = lambda * real(Z) * I using applicators.
%
%   Inputs:
%       eigCount: Number of requested eigenvalues.
%       applicators: Structs of the applicators associated with the problem.
%       options.innerTolerance (default 1e-6): Tolerance of the PCG solver.
%       options.innerMaxIterCount (default 1e3): Maximum number of iterations for the inner PCG solver.
%       options.outerTolerance (default 1e-6): Tolerance passed to eigs function.
%       options.outerMaxIterCount (default 1e3): Maximum number of iterations for the eigs function.
%       options.pcgWarn (default true): Shows warnings if the PCG fails.
%       options.whichEigs (default 'largestabs'): Eigenvalues which to search.
%
%   Outputs:
%       lambda: 1 by k vector of eigenvalues.
%       V: n by k matrix of eigenvector.
%       info: Structure with useful informations about the solution run.

    if nargin < 5
        options = struct();
    end

    innerTolerance = getOption(options, 'innerTolerance', 1e-12);
    innerMaxIterCount = getOption(options, 'innerMaxIterCount', 1e3);
    outerTolerance = getOption(options, 'outerTolerance', 1e-12);
    outerMaxIterCount = getOption(options, 'outerMaxIterCount', 1e3);
    pcgWarn = getOption(options, 'pcgWarn', true);
    whichEigs = getOption(options, 'whichEigs', 'largestabs');

    info.innerCallsCount = 0;
    info.pcgFailCount = 0;
    info.pcgTotalIterCount = 0;

    % Inner solver that solves imag(Z) * I = V for I, we use PCG
    function output = applyInvImagZ(solver, input)
        solverOutput = solver.solve( ...
            'initCurrent', zeros(applicators.imagGreensMatrix.inputSize(), 1), ...
            'excitation', input, ...
            'tolerance', innerTolerance, ...
            'maxIterCount', innerMaxIterCount ...
        );

        output = complex(solverOutput.postRunInfo().current);

        info.pcgTotalIterCount = info.pcgTotalIterCount + solverOutput.count();

        if ~solverOutput.postRunInfo().converged
            if pcgWarn
                WARN('Inner PCG solver does not converged.');
            end

            info.pcgFailCount = info.pcgFailCount + 1;
        end
    end

    % Define the composed operator of the equivalent eigenvalue problem
    function output = applyInvImagZRealZ(solver, input)
        [nrows, ncols] = size(input);

        output = zeros(nrows, ncols, 'like', input);

        for col = 1:ncols
            t = applicators.realGreensMatrix.apply(input(:, col));
            output(:, col) = applyInvImagZ(solver, t);
        end
    end

    % Initialize the solver for the inner part.
    solver = rooftopMoM.core.solver.PreconditionedConjugateGradient(applicators.imagGreensMatrix, applicators.imagPreconditioner);

    % Make a function handler for the eigs function.
    eigsHandle = @(x) applyInvImagZRealZ(solver, x);

    [V, D] = eigs(eigsHandle, applicators.realGreensMatrix.inputSize(), eigCount, whichEigs, struct('tol', outerTolerance, 'maxit', outerMaxIterCount));

    % And finally invert.
    lambda = 1 ./ diag(D);
end
