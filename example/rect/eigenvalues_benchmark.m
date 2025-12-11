%% eigenvalues benchmark
%
%   Description:
%       This benchmark calculates the Characteristic Modes of a patch antenna.
%       It demonstrates the usage of accelerated impedance matrix applicators
%       to reduce computational complexity:
%
%       1. Matrix-Vector Product (Z * I): Reduced from O(N^2) to O(N log N).
%       2. Inverse Operator (Z^(-1) * V): Typically O(N^3), but reduced here 
%           to O(N log N) using iterative solvers and FFT acceleration.
%
%       The script sets up the geometry, prepares the solver context, and solves
%       the Generalized Eigenvalue Problem (GEP) to find the characteristic modes.

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
rooftopCountPerOneWavelength = 500;

% Antennas with electrical size below this threshold will use the same
% discretization as an antenna whose electrical size equals the threshold.
electricalSizeThreshold = 5;

% Antenna side aspect ratio (sizeX:sizeY) will be 1:antennaSidesRatio.
antennaSidesRatio = 2;

% Order of the Gauss–Legendre quadrature.
quadratureOrder = 5;

% Define the incident plane wave.
direction = [0; 0; -1];
polarization = [0; 1; 0];
waveIntensity = 1;

% Memory allocation for a single chunk (number of elements, complex doubles).
chunkSize = 5e6;

% Iterative solver settings.
tolerance = 1e-6;
maxIterCount = 1e4;

%% Benchmark body
[countX, countY, sizeX, sizeY] = generateAntennaDimensions(electricalSize, antennaSidesRatio, rooftopCountPerOneWavelength);

% Compute the wavenumber for the current antenna geometry.
wavenumber = 2 * electricalSize / sqrt(sizeX^2 * countX^2 + sizeY^2 * countY^2);

DEBUG('eigenvalues_benchmark: Preparing shared data for the solvers.');

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

DEBUG('eigenvalues_benchmark: Computing with the FFT preconditioned CG solver');

applicators = generateApplicators(solverContext, ...
    'sizeX', sizeX, ...
    'sizeY', sizeY, ...
    'preconditionerType', 'Strang' ...
);

%% First eigenvalue problem (eig1): imag(Z) * I = lambda * real(Z) * I
[lambda, V, info] = solveEig1Problem(eigCount, ...
    applicators, ...
    struct( ...
        'innerTolerance', tolerance, ...
        'innerMaxIterCount', maxIterCount, ...
        'whichEigs', whichEigs ...
    ) ...
);

%% Display the eigenvalues
disp(lambda);
