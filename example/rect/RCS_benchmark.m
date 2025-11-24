%% RCS benchmark
clc; clearvars;

% Include paths.
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

%% Prepare the benchmark inputs and outputs
% Choose the minimal and maximal electrical size and the sample count
kaMin = 1e-3;
kaMax = 1e-2;
sampleCount = 10;

% Constants for the wavenumber and antenna size calculation.
rooftopCountPerOneWavelength = 10;

% The final antenna ratio will be 1:antennaSidesRatio.
antennaSidesRatio = 2;

% Select the order of the Gauss-Legendre quadrature.
quadratureOrder = 3;

% Define the plane wave
direction = [0; 0; -1];
polarization = [0; 1; 0];
waveIntensity = 1;

% Memory allocation for single chunk in terms of element (complex double) count
chunkSize = 1e7;

% Iterative solver settings
relativeTolerance = 1e-6;
maxIterationCount = 1e4;

output = solveAntennaCurrent('useSplitFFT', false);
output.solverOutput.postRunInfo().current
