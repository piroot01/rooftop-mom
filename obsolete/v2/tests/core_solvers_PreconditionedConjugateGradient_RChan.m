clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

loggerDir = fullfile(currentDir, '..', '3rd', 'matlab-logger', 'src');
addpath(loggerDir);

% Configure the logger.
logger.Logger.configure( ...
    'timeFormat', 'ss.SSS', ...
    'minLevel', logger.LogLevel.TRACE, ...
    'sinks', { ...
        logger.sink.ConsoleSink(), ...
    } ...
);

countX = 100;
countY = 1;
sizeX = 1;
sizeY = 1;

order = 3;
elementCount = 1e6;
k = pi / (sizeY * countY);
V = 1;

quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);
rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);

rooftops = rooftopMoM.basisFunctions.UniformRectPrimitiveUniformRooftops(rectPrimitive, quad);
rooftops.generate();

planner = rooftopMoM.core.greensMatrix.AssemblerChunkPlanner(elementCount);
context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);
assembler = rooftopMoM.core.greensMatrix.FullAssembler(planner);

greensMatrix = assembler.generate(context);

spectralGreensMatrix = greensMatrix.generateSpectralGreensMatrix();

greensMatrixApplicator = rooftopMoM.core.greensMatrix.FFTAcceleratedApplicator( ...
    context, spectralGreensMatrix, greensMatrix.rooftopContent());

% Generate the Strangs preconditioner
preconditioner = rooftopMoM.core.greensMatrix.preconditioner.utils.generatePreconditioner1( ...
    @rooftopMoM.core.greensMatrix.preconditioner.utils.generateTChanPreconditionerKernel, ...
    greensMatrix, sizeX, sizeY);

preconditionerApplicator = rooftopMoM.core.greensMatrix.preconditioner.FFTAcceleratedApplicator( ...
    context, preconditioner, greensMatrix.rooftopContent());
preconditionerApplicator.generate();

solver = rooftopMoM.core.solver.PreconditionedConjugateGradient(greensMatrixApplicator, preconditionerApplicator);

mode = rooftopMoM.operator.utils.generateMiddlePortMode(rooftops.count());
excitation = rooftopMoM.operator.DeltaGapFeederExcitation.generate(rooftops.count(), mode, V);

tolerance = 1e-6;
maxIterCount = 10000;

outputConvSolver = solver.solve( ...
    'excitation', excitation, ...
    'tolerance', tolerance, ...
    'maxIterCount', maxIterCount ...
);

fprintf('Number of the iterations: %d\n', outputConvSolver.count());
