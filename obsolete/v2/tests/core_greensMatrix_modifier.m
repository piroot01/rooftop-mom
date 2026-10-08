clc; clearvars; format long;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

loggerDir = fullfile(currentDir, '..', '3rd', 'matlab-logger', 'src');
addpath(loggerDir);

% Configure the logger.
logger.Logger.configure( ...
    'timeFormat', 'ss.SSS', ...
    'minLevel', logger.LogLevel.DEBUG, ...
    'sinks', { ...
        logger.sink.ConsoleSink(), ...
    } ...
);

rng(0);

countX = 3;
countY = 6;
sizeX = 1;
sizeY = 1;

rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
rect = rooftopMoM.mesh.UniformRect(rectPrimitive);
mesh = rect.generateMesh();

order = 3;

quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);
rooftops = rooftopMoM.basisFunctions.MeshUniformRooftops(mesh, quad);

rooftops.generate();

elementCount = 1e7;

planner = rooftopMoM.core.impedanceMatrix.AssemblerChunkPlanner(elementCount);

k = 1;

context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);

fullAssembler = rooftopMoM.core.impedanceMatrix.FullAssembler(planner);

output = fullAssembler.generateFull(context);

Z = output.value();
V = randi(10, [size(Z, 1), 1]);

% Reference
Iref = imag(Z) \ V

% Iterative solver approach
rooftops = rooftopMoM.basisFunctions.UniformRectPrimitiveUniformRooftops(rectPrimitive, quad);
rooftops.generate();

planner = rooftopMoM.core.greensMatrix.AssemblerChunkPlanner(elementCount);
context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);
assembler = rooftopMoM.core.greensMatrix.FullAssembler(planner);

greensMatrix = assembler.generate(context);

spectralGreensMatrix = greensMatrix.generateSpectralGreensMatrix( ...
    struct('modifier', @imag) ...
);

greensMatrixApplicator = rooftopMoM.core.greensMatrix.FFTAcceleratedApplicator( ...
    context, spectralGreensMatrix, greensMatrix.rooftopContent());

preconditioner = rooftopMoM.core.greensMatrix.preconditioner.utils.generatePreconditioner1( ...
    @rooftopMoM.core.greensMatrix.preconditioner.utils.generateTChanPreconditionerKernel, ...
    greensMatrix, sizeX, sizeY, struct('modifier', @imag));

preconditionerApplicator = rooftopMoM.core.greensMatrix.preconditioner.FFTAcceleratedApplicator( ...
    context, preconditioner, greensMatrix.rooftopContent());
preconditionerApplicator.generate();

solver = rooftopMoM.core.solver.PreconditionedConjugateGradient(greensMatrixApplicator, preconditionerApplicator);

tolerance = 1e-12;
maxIterCount = 1000;

outputConvSolver = solver.solve( ...
    'initCurrent', zeros(rooftops.count(), 1), ...
    'excitation', V, ...
    'tolerance', tolerance, ...
    'maxIterCount', maxIterCount ...
);

fprintf('Number of the iterations: %d\n', outputConvSolver.count());

currentConv = outputConvSolver.postRunInfo().current;

norm(Iref - currentConv)
