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
% The mesh is specified by the number of faces in each direction.
countX = 500;
countY = 1;

% Size of each elementary cell (face) in the mesh.
sizeX = 1;
sizeY = 1;

%% MoM setup
% Order of the Gauss–Legendre quadrature.
% Distance matrix batch size (as a double count).
order = 3;

% Distance matrix batch size in as a double count.
elementCount = 1e7;

% Wavelength of the excitation field (half-wavelength dipole),
% hence k = 2 * pi / lambda = pi / (sizeX * countX).
k = pi / (sizeX * countX);

% Voltage (Volts) of the delta-gap excitation.
V = 1;

% Iterative solver settings: tolerance and maximum number of iterations.
tolerance = 1e-9;
maxIterCount = 1000;

% Create the excitation vector for the delta-gap feeder.
count = rooftopMoM.basisFunctions.utils.generateRooftopCount(countX, countY);
mode = rooftopMoM.operator.utils.generateMiddlePortMode(count);
excitation = rooftopMoM.operator.DeltaGapFeederExcitation.generate(count, mode, V);

%% Naive MoM implementation:
% Generate a primitive description of the mesh. This is a lightweight
% structure holding only the basic tiling info (no vertices yet).
rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);

% Convert the primitive into a 3D-oriented rectangle object.
rect = rooftopMoM.mesh.UniformRect(rectPrimitive);

% Generate the actual mesh of the antenna.
mesh = rect.generateMesh();

% Quadrature for two rectangular faces (Gauss-Legendre).
quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);

% Build rooftop basis on the mesh (also computes any auxiliary info needed).
rooftops = rooftopMoM.basisFunctions.MeshUniformRooftops(mesh, quad);
rooftops.generate();

% Planner controls how chunks are arranged during impedance assembly.
planner = rooftopMoM.core.impedanceMatrix.AssemblerChunkPlanner(elementCount);

% Bundle required data for the assembler.
context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);

% FullAssembler efficiently composes magnetic and electric parts for us.
fullAssembler = rooftopMoM.core.impedanceMatrix.FullAssembler(planner);

strategyNaive = tic;

% Assemble the impedance matrix.
outputNaive = fullAssembler.generateFull(context);

% Using the build-in MATLAB solver find the currents.
currentNaive = rooftopMoM.core.solver.Inversion.solve(outputNaive.value(), excitation);

timeNaive = toc(strategyNaive);

%% Interaction-symmetries-accelerated MoM implementation
% Generate a primitive description of the mesh. This is a lightweight
% structure holding only the basic tiling info (no vertices yet).
rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);

% Convert the primitive into a 3D-oriented rectangle object.
rect = rooftopMoM.mesh.UniformRect(rectPrimitive);

% Generate the actual mesh of the antenna.
mesh = rect.generateMesh();

% Quadrature for two rectangular faces (Gauss-Legendre).
quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);

% Build rooftop basis on the mesh (also computes any auxiliary info needed).
rooftops = rooftopMoM.basisFunctions.MeshUniformRooftops(mesh, quad);
rooftops.generate();

% Planner controls how chunks are arranged during impedance assembly.
planner = rooftopMoM.core.impedanceMatrix.AssemblerChunkPlanner(elementCount);

% Bundle required data for the assembler.
context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);

% FullAssembler efficiently composes magnetic and electric parts for us.
fullAssembler = rooftopMoM.core.impedanceMatrix.FullAssembler(planner);

% InteractionSymmetriesAccelerator wraps the assembler and exploits symmetry
% in the rooftop interactions. This can reduce computation and memory by
% reusing symmetric blocks instead of computing them repeatedly.
accelerator = rooftopMoM.core.impedanceMatrix.InteractionSymmetriesAccelerator(fullAssembler);

strategyAccelerated = tic;

% Assemble the impedance matrix.
outputAccelerated = accelerator.generateFull(context);

% Using the build-in MATLAB solver find the currents.
currentAccelerated = rooftopMoM.core.solver.Inversion.solve(outputAccelerated.value(), excitation);

timeAccelerated = toc(strategyAccelerated);

%% Green's operator via convolution MoM implementation - using conjugate gradient method
% Generate a primitive description of the mesh. This is a lightweight
% structure holding only the basic tiling info (no vertices yet).
rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);

% Quadrature for two rectangular faces (Gauss-Legendre).
quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);

% Build rooftop basis on the mesh (also computes any auxiliary info needed).
rooftops = rooftopMoM.basisFunctions.UniformRectPrimitiveUniformRooftops(rectPrimitive, quad);
rooftops.generate();

% Planner controls how chunks are arranged during impedance assembly.
planner = rooftopMoM.core.greensMatrix.AssemblerChunkPlanner(elementCount);

% Bundle required data for the assembler.
context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);

% FullAssembler efficiently composes magnetic and electric parts for us.
fullAssembler = rooftopMoM.core.greensMatrix.FullAssembler(planner);

strategyConjugateGradient = tic;

% Assemble the greens matrix kernels and prepare the split-FFT applicator.
greensMatrix = fullAssembler.generate(context);

greensMatrixApplicator = rooftopMoM.core.greensMatrix.SplitFFTAcceleratedApplicator( ...
    context, greensMatrix.generateSplitSpectralGreensMatrix(struct('faceSizeX', sizeX, 'faceSizeY', sizeY)), greensMatrix.rooftopContent());
greensMatrixApplicator.generate();

conjugateGradientSolver = rooftopMoM.core.solver.ConjugateGradient(greensMatrixApplicator);

% Run the solver.
% The current and other info is stored within the outputSolver.
outputSolver = conjugateGradientSolver.solve( ...
    'excitation', excitation, ...
    'tolerance', tolerance, ...
    'maxIterCount', maxIterCount ...
);

currentConjugateGradient = outputSolver.postRunInfo().current;

timeConjugateGradient = toc(strategyConjugateGradient);

%% Green's operator via convolution MoM implementation - using preconditioned conjugate gradient method
% Generate a primitive description of the mesh. This is a lightweight
% structure holding only the basic tiling info (no vertices yet).
rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);

% Quadrature for two rectangular faces (Gauss-Legendre).
quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);

% Build rooftop basis on the mesh (also computes any auxiliary info needed).
rooftops = rooftopMoM.basisFunctions.UniformRectPrimitiveUniformRooftops(rectPrimitive, quad);
rooftops.generate();

% Planner controls how chunks are arranged during impedance assembly.
planner = rooftopMoM.core.greensMatrix.AssemblerChunkPlanner(elementCount);

% Bundle required data for the assembler.
context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);

% FullAssembler efficiently composes magnetic and electric parts for us.
fullAssembler = rooftopMoM.core.greensMatrix.FullAssembler(planner);

strategyPreconditionedConjugateGradient = tic;

% Assemble the greens matrix kernels and prepare the split-FFT applicator.
greensMatrix = fullAssembler.generate(context);

greensMatrixApplicator = rooftopMoM.core.greensMatrix.SplitFFTAcceleratedApplicator( ...
    context, greensMatrix.generateSplitSpectralGreensMatrix(struct('faceSizeX', sizeX, 'faceSizeY', sizeY)), greensMatrix.rooftopContent());
greensMatrixApplicator.generate();

% The Strang's circulant approximation of the greens matrix is inverted in the spectral domain.
preconditioner = rooftopMoM.core.greensMatrix.preconditioner.utils.generatePreconditioner1( ...
    @rooftopMoM.core.greensMatrix.preconditioner.utils.generateStrangPreconditionerKernel, ...
    greensMatrix, sizeX, sizeY);

preconditionerApplicator = rooftopMoM.core.greensMatrix.preconditioner.FFTAcceleratedApplicator( ...
    context, preconditioner, greensMatrix.rooftopContent());
preconditionerApplicator.generate();

spectrallyPreconditionedConjugateGradientSolver = rooftopMoM.core.solver.PreconditionedConjugateGradient(greensMatrixApplicator, preconditionerApplicator);

% Run the solver.
% The current and other info is stored within the outputSolver.
outputSolver = spectrallyPreconditionedConjugateGradientSolver.solve( ...
    'excitation', excitation, ...
    'tolerance', tolerance, ...
    'maxIterCount', maxIterCount ...
);

currentPreconditionedConjugateGradient = outputSolver.postRunInfo().current;

timePreconditionedConjugateGradient = toc(strategyPreconditionedConjugateGradient);

%% Summary table (append at end of script)
% Prints a compact Linux-style ASCII table summarizing each MoM strategy.

% Prepare rows
names = {'Naive', 'Accelerated', 'Conjugate CG', 'Precond. Conj. CG'};
methods = { ...
    'Full MoM assembly', ...
    'Symmetry-accelerated assembly', ...
    'Green''s operator via convolution', ...
    'Green''s operator + spectral preconditioner' ...
};
solvers = { ...
    'Built-in direct solver', ...
    'Built-in direct solver', ...
    'Conjugate Gradient (iterative)', ...
    'Spectrally preconditioned CG (iterative)' ...
};

% Collect times.
times = [ ...
    timeNaive, ...
    timeAccelerated, ...
    timeConjugateGradient, ...
    timePreconditionedConjugateGradient ...
];

% Print ASCII table
line = '+----------------------+-----------------------------------------------+------------------------------------------+------------+';
fprintf('%s\n', line);
fprintf('| %-20s | %-45s | %-40s | %-10s |\n', 'Strategy', 'Method', 'Solver', 'Time (s)');
fprintf('%s\n', line);
for i = 1:length(names)
    t = times(i);
    if isnan(t)
        tstr = '   N/A    ';
    else
        tstr = sprintf('%10.4f', t);
    end
    fprintf('| %-20s | %-45s | %-40s | %s |\n', names{i}, methods{i}, solvers{i}, tstr);
end
fprintf('%s\n', line);
