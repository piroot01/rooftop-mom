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
    'minLevel', logger.LogLevel.DEBUG, ...
    'sinks', { ...
        logger.sink.ConsoleSink(), ...
    } ...
);

countX = 100;
countY = 100;
sizeX = 0.5;
sizeY = 1.7;

rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
rect = rooftopMoM.mesh.UniformRect(rectPrimitive);
mesh = rect.generateMesh();

order = 3;
elementCount = 1e7;
k = pi / (sizeY * countY);
V = 1;

quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);

%
% Classical naive approach
%

rooftops = rooftopMoM.basisFunctions.MeshUniformRooftops(mesh, quad);
rooftops.generate();

%planner = rooftopMoM.core.impedanceMatrix.AssemblerChunkPlanner(elementCount);
%context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);
%fullAssembler = rooftopMoM.core.impedanceMatrix.FullAssembler(planner);

%outputNaive = fullAssembler.generateFull(context);

mode = rooftopMoM.operators.utils.generateMiddlePortMode(rooftops.count());
excitation = rooftopMoM.operators.DeltaGapFeederExcitation.generate(rooftops.count(), mode, V);

%currentNaive = rooftopMoM.core.solver.Inversion.solve(outputNaive.value(), excitation);

%
% Convolutional approach
%

rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
rooftops = rooftopMoM.basisFunctions.UniformRectPrimitiveUniformRooftops(rectPrimitive, quad);
rooftops.generate();

planner = rooftopMoM.core.greensMatrix.AssemblerChunkPlanner(elementCount);
context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);

assembler = rooftopMoM.core.greensMatrix.FullAssembler(planner);

greensMatrix = assembler.generateFull(context);

greensMatrixApplicator = rooftopMoM.core.greensMatrix.FFTAcceleratedApplicator(context, greensMatrix);
greensMatrixApplicator.generate();

conjugateGradientSolver = rooftopMoM.core.solver.ConjugateGradient(greensMatrixApplicator);

tolerance = 1e-3;
maxIterCount = 2000;

outputConvSolver = conjugateGradientSolver.solve( ...
    'excitation', excitation, ...
    'tolerance', tolerance, ...
    'maxIterCount', maxIterCount ...
);

fprintf('Number of the iterations: %d\n', outputConvSolver.count());
currentConv = outputConvSolver.postRunInfo().current;

greensMatrixApplicator = rooftopMoM.core.greensMatrix.SplitFFTAcceleratedApplicator(context, greensMatrix);
greensMatrixApplicator.generate();

conjugateGradientSolver = rooftopMoM.core.solver.ConjugateGradient(greensMatrixApplicator);

tolerance = 1e-3;
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