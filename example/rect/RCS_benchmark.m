%% RCS benchmark
%
%   Description:
%       This benchmark compares accuracy and solving time of two iterative solvers
%       that use the conjugate-gradient method with a spectral preconditioner but
%       different Green's-matrix applicators. The FFTMult variant performs the
%       matrix–vector product using classical FFT-accelerated convolution; the
%       splitFFTMult variant partitions arrays to avoid zero-padding, reducing memory
%       usage and becoming faster for large problems. We evaluate both methods by
%       computing the radar cross section (RCS) of a rectangular antenna (side ratio
%       1:2) illuminated by a plane wave with k-vector normal to the antenna and
%       polarization along the longer edge (unit amplitude). The first plot shows
%       RCS normalized by antenna area as a function of electrical size; the second
%       plot shows solver run times versus electrical size.

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

%% Enable caching of rooftop data arrays to speed up operator evaluations
rooftopMoM.Cache.setStatus(true);
rooftopMoM.Cache.setCacheAll(true);

%% Prepare benchmark inputs and outputs
% Choose minimum and maximum electrical sizes and the number of samples.
kaMin = 1e-1;
kaMax = 100;
sampleCount = 20;

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
chunkSize = 1e7;

% Iterative solver settings.
tolerance = 1e-6;
maxIterCount = 1e4;

%% Main benchmark loop
ka = logspace(log10(kaMin), log10(kaMax), sampleCount);

% Preallocate a cell array to store per-sample results (each entry is a struct).
benchmarkOutput = cell(1, sampleCount);

for n = 1:sampleCount
    DEBUG('RCS_benchmark: Executing iteration %g out of %g.', n, sampleCount);
    electricalSize = ka(n);

    % Choose discretization: use threshold if the current electrical size is small.
    if electricalSize < electricalSizeThreshold 
        [countX, countY, sizeX, sizeY] = generateAntennaDimensions(electricalSizeThreshold, antennaSidesRatio, rooftopCountPerOneWavelength);
    else
        [countX, countY, sizeX, sizeY] = generateAntennaDimensions(electricalSize, antennaSidesRatio, rooftopCountPerOneWavelength);
    end

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

    DEBUG('RCS_benchmark: Computing with the FFT-preconditioned CG solver');

    result_FFTMult = solveAntennaCurrent(solverContext, ...
        'tolerance', tolerance, ...
        'maxIterCount', maxIterCount, ...
        'useSplitFFT', false ...
    );

    current = result_FFTMult.solverOutput.postRunInfo().current;

    rho_FFTMult = computeRCS(wavenumber, current, solverContext.excitation, waveIntensity);
    clear current;

    DEBUG('RCS_benchmark: Computing with the split-FFT preconditioned CG solver');

    result_splitFFTMult = solveAntennaCurrent(solverContext, ...
        'tolerance', tolerance, ...
        'maxIterCount', maxIterCount, ...
        'useSplitFFT', true ...
    );

    current = result_splitFFTMult.solverOutput.postRunInfo().current;

    rho_splitFFTMult = computeRCS(wavenumber, current, solverContext.excitation, waveIntensity);
    clear current;

    % Retrieve runtime information.
    solvingTime_FFTMult = result_FFTMult.solverOutput.postRunInfo().solvingTime;
    solvingTime_splitFFTMult = result_splitFFTMult.solverOutput.postRunInfo().solvingTime;
    DEBUG('RCS_benchmark: FFTMult variant took %g s and splitFFTMult took %g s (%gx speed up)', solvingTime_FFTMult, solvingTime_splitFFTMult, solvingTime_FFTMult / solvingTime_splitFFTMult);

    % Store results for this sample.
    benchmarkOutput{n} = struct( ...
        'wavenumber', wavenumber, ...
        'antennaArea', countX * countY, ...
        'RCS_FFTMult', rho_FFTMult, ...
        'RCS_splitFFTMult', rho_splitFFTMult, ...
        'solvingTime_FFTMult', solvingTime_FFTMult, ...
        'solvingTime_splitFFTMult', solvingTime_splitFFTMult ...
    );
end

%% Data analysis
% Extract the benchmark results into vectors.
wavenumber = cellfun(@(c) c.wavenumber, benchmarkOutput);
antennaArea = cellfun(@(c) c.antennaArea, benchmarkOutput);

RCS_FFTMult = cellfun(@(c) c.RCS_FFTMult, benchmarkOutput);
RCS_splitFFTMult = cellfun(@(c) c.RCS_splitFFTMult, benchmarkOutput);
RCS_norm_FFTMult = RCS_FFTMult ./ antennaArea;
RCS_norm_splitFFTMult = RCS_splitFFTMult ./ antennaArea;

solvingTime_FFTMult = cellfun(@(c) c.solvingTime_FFTMult, benchmarkOutput);
solvingTime_splitFFTMult = cellfun(@(c) c.solvingTime_splitFFTMult, benchmarkOutput);

% Compare normalized RCS for both solvers.
figure('Name', 'RCS_benchmark: normalized RCS comparison');
loglog(ka, RCS_norm_FFTMult, 'x--', 'LineWidth', 2, 'MarkerSize', 10);
hold on;
loglog(ka, RCS_norm_splitFFTMult, 'd', 'LineWidth', 2, 'MarkerSize', 10);
loglog(ka, wavenumber.^2 .* antennaArea / pi, '--', 'LineWidth', 2);
hold off;

xlim([kaMin, kaMax]);
xlabel('ka = 2\pi / \lambda^2 [-]')
ylabel('\sigma / A [-]');
legend({'FFTMult', 'splitFFTMult', '4\pi / \lambda^2'}, 'Location', 'best');
grid on;

% Compare solver runtimes.
figure('Name', 'RCS_benchmark: solving time comparison');
loglog(ka, solvingTime_FFTMult, 'x--', 'LineWidth', 2, 'MarkerSize', 10);
hold on;
loglog(ka, solvingTime_splitFFTMult, 'x--', 'LineWidth', 2, 'MarkerSize', 10);
hold off;

xlim([kaMin, kaMax]);
xlabel('ka = 2\pi / \lambda^2 [-]')
ylabel('t [s]');
legend({'FFTMult', 'splitFFTMult'}, 'Location', 'best');
grid on;
