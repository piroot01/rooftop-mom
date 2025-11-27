clc; clearvars;

%profile clear;
%profile -memory on;

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

countX = 200;
countY = 100;
sizeX = 1.5;
sizeY = 0.5;

order = 2;
elementCount = 10000;
k = pi / (sizeY * countY);
V = 1;

quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);

%
% Convolutional approach
%

rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
rooftops = rooftopMoM.basisFunctions.UniformRectPrimitiveUniformRooftops(rectPrimitive, quad);

tic;
rooftops.generate();
toc;

rooftops.count()

planner = rooftopMoM.core.greensMatrix.AssemblerChunkPlanner(elementCount);
context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);

assembler = rooftopMoM.core.greensMatrix.FullAssembler(planner);

tic;
greensMatrix = assembler.generateFull(context);
toc;

greensMatrixApplicator = rooftopMoM.core.greensMatrix.FFTAcceleratedApplicator(context, greensMatrix);
greensMatrixApplicator.generate();

conjugateGradientSolver = rooftopMoM.core.solver.ConjugateGradient(greensMatrixApplicator);

mode = rooftopMoM.operator.utils.generateMiddlePortMode(rooftops.count());
excitation = rooftopMoM.operator.DeltaGapFeederExcitation.generate(rooftops.count(), mode, V);

tolerance = 1e-1;
maxIterCount = 10000;

outputConvSolver = conjugateGradientSolver.solve( ...
    'excitation', excitation, ...
    'tolerance', tolerance, ...
    'maxIterCount', maxIterCount ...
);

fprintf('Number of the iterations: %d\n', outputConvSolver.count());
currentConv = outputConvSolver.postRunInfo().current;

%norm(currentConv - currentNaive)

greensMatrixApplicator = rooftopMoM.core.greensMatrix.SplitFFTAcceleratedApplicator(context, greensMatrix);
greensMatrixApplicator.generate();

conjugateGradientSolver = rooftopMoM.core.solver.ConjugateGradient(greensMatrixApplicator);

tolerance = 1e-1;
maxIterCount = 2000;

outputConvSolver = conjugateGradientSolver.solve( ...
    'excitation', excitation, ...
    'tolerance', tolerance, ...
    'maxIterCount', maxIterCount ...
);

fprintf('Number of the iterations: %d\n', outputConvSolver.count());
currentConvSplit = outputConvSolver.postRunInfo().current;

% check the norm
norm(currentConvSplit - currentConv)