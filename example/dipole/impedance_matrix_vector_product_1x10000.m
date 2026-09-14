clc; clearvars;

% Include paths.
currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', '..', 'src');
addpath(srcDir);

loggerDir = fullfile(currentDir, '..', '..', '3rd', 'matlab-logger', 'src');
addpath(loggerDir);

% Configure the logger.
logger.Logger.configure( ...
    'timeFormat', 'ss.SSS', ...
    'minLevel', logger.LogLevel.INFO, ...
    'sinks', { ...
        logger.sink.ConsoleSink(), ...
    } ...
);

%% Antenna mesh setup:
% The mesh is specified by the number of faces in each direction.
countX = 10000;
countY = 1;

% Size of each elementary cell (face) in the mesh.
sizeX = 1;
sizeY = 1;

% Generate a primitive description of the mesh. This is a lightweight
% structure holding only the basic tiling info (no vertices yet).
rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);

%% MoM setup:
% Order of the Gauss–Legendre quadrature.
order = 3;

% Distance matrix batch size in as a double count.
elementCount = 1e6;

% Wavelength of the excitation field (half-wavelength dipole),
% hence k = 2 * pi / lambda = pi / (sizeX * countX).
k = pi / (sizeX * countX);

% Number of the repeated applications used to measure the mean application time.
repetitionCount = 100;

% Quadrature for two rectangular faces (Gauss-Legendre).
quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);

%% Raw impedance matrix:
% Convert the primitive into a 3D-oriented rectangle object.
rect = rooftopMoM.mesh.UniformRect(rectPrimitive);

% Generate the actual mesh of the antenna.
mesh = rect.generateMesh();

% Build rooftop basis on the mesh (also computes any auxiliary info needed).
meshRooftops = rooftopMoM.basisFunctions.MeshUniformRooftops(mesh, quad);
meshRooftops.generate();

% Planner controls how chunks are arranged during impedance assembly.
planner = rooftopMoM.core.impedanceMatrix.AssemblerChunkPlanner(elementCount);

% Bundle required data for the assembler.
context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, meshRooftops, k);

% InteractionSymmetriesAccelerator wraps the FullAssembler and reuses the
% symmetric blocks of the rooftop interactions.
fullAssembler = rooftopMoM.core.impedanceMatrix.FullAssembler(planner);
accelerator = rooftopMoM.core.impedanceMatrix.InteractionSymmetriesAccelerator(fullAssembler);

tic;
impedanceMatrix = accelerator.generateFull(context).value();
elapsed = toc;
fprintf('Impedance matrix assembly finished in %.2f s.\n', elapsed);

%% Compressed greens matrix:
% Build rooftop basis on the primitive (no mesh is needed).
rooftops = rooftopMoM.basisFunctions.UniformRectPrimitiveUniformRooftops(rectPrimitive, quad);
rooftops.generate();

% Planner controls how chunks are arranged during greens matrix assembly.
greensPlanner = rooftopMoM.core.greensMatrix.AssemblerChunkPlanner(elementCount);

% Bundle required data for the assembler.
greensContext = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);

% FullAssembler stores only the Toeplitz kernels, the full matrix is never formed.
greensAssembler = rooftopMoM.core.greensMatrix.FullAssembler(greensPlanner);

tic;
greensMatrix = greensAssembler.generate(greensContext);
elapsed = toc;
fprintf('Greens matrix assembly finished in %.2f s.\n', elapsed);

% The face sizes enter the spectra through the electric part.
spectralOptions = struct('faceSizeX', sizeX, 'faceSizeY', sizeY);

% FFT accelerated applicator (kernels embedded into zero-padded circulants).
fftApplicator = rooftopMoM.core.greensMatrix.FFTAcceleratedApplicator( ...
    greensContext, greensMatrix.generateSpectralGreensMatrix(spectralOptions), greensMatrix.rooftopContent());
fftApplicator.generate();

% Split-FFT accelerated applicator (subband spectra of the grid size).
splitFFTApplicator = rooftopMoM.core.greensMatrix.SplitFFTAcceleratedApplicator( ...
    greensContext, greensMatrix.generateSplitSpectralGreensMatrix(spectralOptions), greensMatrix.rooftopContent());
splitFFTApplicator.generate();

%% Random input vector:
basisSize = rooftops.count();
randomVector = randn(basisSize, 1) + 1i * randn(basisSize, 1);

%% Raw impedance matrix application:
% Warm-up application.
outputRaw = impedanceMatrix * randomVector;

tic;
for n = 1:repetitionCount
    outputRaw = impedanceMatrix * randomVector;
end
timeRaw = toc / repetitionCount;

%% FFT accelerated application:
% Warm-up application.
outputFFT = fftApplicator.apply(randomVector);

tic;
for n = 1:repetitionCount
    outputFFT = fftApplicator.apply(randomVector);
end
timeFFT = toc / repetitionCount;

%% Split-FFT accelerated application:
% Warm-up application.
outputSplitFFT = splitFFTApplicator.apply(randomVector);

tic;
for n = 1:repetitionCount
    outputSplitFFT = splitFFTApplicator.apply(randomVector);
end
timeSplitFFT = toc / repetitionCount;

%% Summary table
% Prints a compact Linux-style ASCII table summarizing each application strategy.

% Prepare rows
names = {'Raw', 'FFT', 'Split FFT'};
methods = { ...
    'Dense impedance matrix product', ...
    'Zero-padded circulant embedding', ...
    'Subband spectra of the grid size' ...
};

% Collect times and relative errors with respect to the raw product.
times = [timeRaw, timeFFT, timeSplitFFT];
errors = [ ...
    0, ...
    norm(outputFFT - outputRaw) / norm(outputRaw), ...
    norm(outputSplitFFT - outputRaw) / norm(outputRaw) ...
];

% Print ASCII table
fprintf('Applying the impedance matrix of %d rooftops.\n', basisSize);

line = '+------------+----------------------------------+------------+------------+------------+';
fprintf('%s\n', line);
fprintf('| %-10s | %-32s | %-10s | %-10s | %-10s |\n', 'Strategy', 'Method', 'Time (ms)', 'Speed-up', 'Rel. error');
fprintf('%s\n', line);
for i = 1:length(names)
    fprintf('| %-10s | %-32s | %10.4f | %10.2f | %10.2e |\n', names{i}, methods{i}, 1e3 * times(i), times(1) / times(i), errors(i));
end
fprintf('%s\n', line);
