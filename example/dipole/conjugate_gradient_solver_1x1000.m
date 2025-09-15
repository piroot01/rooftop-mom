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
    'minLevel', logger.LogLevel.DEBUG, ...
    'sinks', { ...
        logger.sink.ConsoleSink(), ...
    } ...
);

%% Antenna mesh setup:
% The mesh is described by how many faces are in each direction.
countX = 1000;
countY = 1;

% Size of each elementary cell (face) of the mesh.
sizeX = 1;
sizeY = 1;

% Generate a primitive description of the mesh. This is a lightweight
% structure holding only the basic tiling info (no vertices yet).
rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);

%% Quadrature / rooftop basis setup:
% Order controls accuracy of Gauss-Legendre quadrature. Higher => more accurate,
% but memory and CPU cost grow quickly (roughly O(N^4) with order N).
order = 5;

% Quadrature for two rectangular faces (Gauss-Legendre).
quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);

% Build rooftop basis on the mesh (also computes any auxiliary info needed).
rooftops = rooftopMoM.basisFunctions.UniformRectPrimitiveUniformRooftops(rectPrimitive, quad);
rooftops.generate();

%% Assembly planning and wavenumber:
% How many elements of the distance matrix to compute per chunk.
elementCount = 1e6;

% Planner controls how chunks are arranged during impedance assembly.
planner = rooftopMoM.core.greensMatrix.AssemblerChunkPlanner(elementCount);

% Example wavenumber: lambda = 2 * sizeX * countX (lambda-half dipole),
% so k = 2 * pi / lambda = pi / (sizeX * countX).
k = pi / (sizeX * countX);

% Bundle required data for the assembler.
context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);

%% Assemble full greens matrix:
% FullAssembler efficiently composes magnetic and electric parts for us.
fullAssembler = rooftopMoM.core.greensMatrix.FullAssembler(planner);

% Start assembly (result stored in 'output'; matrix accessible via output.value())
tic;
output = fullAssembler.generateFull(context);
elapsed = toc;
fprintf('Assembly finished in %.2f s.\n', elapsed);

%% Currents for delta-gap feeder.
% Compute the surface-current distribution produced by a delta-gap feed.
% We use a unit-voltage excitation (V = 1 V) applied at the center port and
% solve for rooftop-basis coefficients (these coefficients represent the
% surface currents).

% Find the currents distributed on the antenna surface.
% Create the excitation vector of the delta-gap feeder with an input voltage V = 1 V.
V = 1;

mode = rooftopMoM.operators.utils.generateMiddlePortMode(rooftops);
excitation = rooftopMoM.operators.DeltaGapFeederExcitation.generate(rooftops, mode, V);

conjugateGradientSolver = rooftopMoM.core.solver.ConjugateGradient(context, output);

% We choose the tolerance and the maximal number of iterations that the solver does.
tolerance = 1e-4;
maxIterCount = 2000;

% Run the solver.
% The current and other info is stored within the outputSolver.
outputSolver = conjugateGradientSolver.solve( ...
    'excitation', excitation, ...
    'tolerance', tolerance, ...
    'maxIterCount', maxIterCount ...
);

fprintf('Solving finished in %.2f s.\n', outputSolver.postRunInfo().solvingTime);
