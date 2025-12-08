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
countY = 400;
sizeX = 1;
sizeY = 1;

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
%{
rooftops = rooftopMoM.basisFunctions.MeshUniformRooftops(mesh, quad);
rooftops.generate();

planner = rooftopMoM.core.impedanceMatrix.AssemblerChunkPlanner(elementCount);
context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);
fullAssembler = rooftopMoM.core.impedanceMatrix.FullAssembler(planner);

outputNaive = fullAssembler.generateFull(context);

mode = rooftopMoM.operator.utils.generateMiddlePortMode(rooftops.count());
excitation = rooftopMoM.operator.DeltaGapFeederExcitation.generate(rooftops.count(), mode, V);

currentNaive = rooftopMoM.core.solver.Inversion.solve(outputNaive.value(), excitation);
%}
%
% Convolutional approach
%

rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
rooftops = rooftopMoM.basisFunctions.UniformRectPrimitiveUniformRooftops(rectPrimitive, quad);
rooftops.generate();

mode = rooftopMoM.operator.utils.generateMiddlePortMode(rooftops.count());
excitation = rooftopMoM.operator.DeltaGapFeederExcitation.generate(rooftops.count(), mode, V);

planner = rooftopMoM.core.greensMatrix.AssemblerChunkPlanner(elementCount);
context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);

assembler = rooftopMoM.core.greensMatrix.FullAssembler(planner);

greensMatrix = assembler.generateFull(context);

spectralGreensMatrix = rooftopMoM.core.greensMatrix.SpectralMatrixData(context);
spectralGreensMatrix.generate(greensMatrix);

% Greens matrix applicator: Z * I = V
greensMatrixApplicator = rooftopMoM.core.greensMatrix.FFTAcceleratedApplicator(context, spectralGreensMatrix);
greensMatrixApplicator.generate();

strangPreconditioner = rooftopMoM.core.greensMatrix.preconditioner.StrangMatrixData(context);
strangPreconditioner.generate(greensMatrix);

preconditionerApplicator = rooftopMoM.core.greensMatrix.preconditioner.FFTAcceleratedCirculant(context, strangPreconditioner);

conjugateGradientSolver = rooftopMoM.core.solver.PreconditionedConjugateGradient(greensMatrixApplicator, preconditionerApplicator);

tolerance = 1e-6;
maxIterCount = 10;

outputConvSolver = conjugateGradientSolver.solve( ...
    'excitation', excitation, ...
    'tolerance', tolerance, ...
    'maxIterCount', maxIterCount ...
);

fprintf('Number of the iterations: %d\n', outputConvSolver.count());
currentConv = outputConvSolver.postRunInfo().current;

preconditionerApplicator = rooftopMoM.core.greensMatrix.preconditioner.FFTAcceleratedToeplitz(context, spectralGreensMatrix);
preconditionerApplicator.generate();

conjugateGradientSolver = rooftopMoM.core.solver.PreconditionedConjugateGradient(greensMatrixApplicator, preconditionerApplicator);

tolerance = 1e-6;
maxIterCount = 2000;

outputConvSolver = conjugateGradientSolver.solve( ...
    'excitation', excitation, ...
    'tolerance', tolerance, ...
    'maxIterCount', maxIterCount ...
);

fprintf('Number of the iterations: %d\n', outputConvSolver.count());
currentConv2 = outputConvSolver.postRunInfo().current;

%
% Compare the currents.
%

%norm(currentConv - currentNaive)
%norm(currentConv2 - currentNaive)
norm(currentConv - currentConv2)
