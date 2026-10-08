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
    'minLevel', logger.LogLevel.DEBUG, ...
    'sinks', { ...
        logger.sink.ConsoleSink(), ...
    } ...
);

%% Prepare benchmark inputs
electricalSize = 10;

% Number of requested eigenvalues
eigCount = 10;
whichEigs = 'largestabs';

% Constants used to compute the wavenumber and antenna size.
% Number of rooftops per wavelength to aim for (10 keeps the problem below 1000 rooftops).
rooftopCountPerOneWavelength = 10;

% Antennas with electrical size below this threshold will use the same
% discretization as an antenna whose electrical size equals the threshold.
electricalSizeThreshold = 5;

% Antenna side aspect ratio (sizeX:sizeY) will be 1:antennaSidesRatio.
antennaSidesRatio = 2;

% Order of the Gauss–Legendre quadrature.
quadratureOrder = 10;

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

%% Comparison with reference data
lambdaRef = [ ...
   0.003674276647855e2; ...
  -0.020072856748803e2; ...
   0.053199894297259e2; ...
  -0.071408911064240e2; ...
  -0.089574441948060e2; ...
   0.304952611989058e2; ...
  -1.686742322387722e2; ...
  -1.796098111244691e2; ...
   2.119329913052550e2; ...
  -2.923984125245223e2; ...
];

toLogScale = @(x) 10*log10(abs(x));

figure
title(sprintf('Characteristic modes accuracy comparison - %d rooftops', applicators.realGreensMatrix.inputSize()));
hold on;
plot(toLogScale(lambda), '--d', 'LineWidth', 2, 'MarkerSize', 8);
plot(toLogScale(lambdaRef), '--x', 'LineWidth', 2, 'MarkerSize', 8);
hold off;
xlabel('modal index n [-]');
ylabel('10log_{10}(|\lambda_n|) [-]')
legend('\lambda_n', 'reference \lambda_n');
