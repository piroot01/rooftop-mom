clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
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

countX = 1;
countY = 20;
sizeX = 1;
sizeY = 1;

rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
rect = rooftopMoM.mesh.UniformRect(rectPrimitive);
mesh = rect.generateMesh();

order = 5;
elementCount = 100000;
k = pi / (sizeY * countY);
V = 1;

quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);

%
% Classical naive approach
%

rooftopCount = rooftopMoM.basisFunctions.utils.generateRooftopCount(countX, countY);
genome = true(rooftopCount, 1);
genome([1, end]) = false

rooftops = rooftopMoM.basisFunctions.MeshUniformRooftopsGenome(mesh, quad, genome);
rooftops.generate();

planner = rooftopMoM.core.impedanceMatrix.AssemblerChunkPlanner(elementCount);
context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);
fullAssembler = rooftopMoM.core.impedanceMatrix.FullAssembler(planner);

outputNaive = fullAssembler.generateFull(context);

mode = rooftopMoM.operators.utils.generateMiddlePortMode(rooftops);
excitation = rooftopMoM.operators.DeltaGapFeederExcitation.generate(rooftops, mode, V);

currentNaive = rooftopMoM.core.solver.Inversion.solve(outputNaive.value(), excitation)

%
% Convolutional approach
%

rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
rooftops = rooftopMoM.basisFunctions.UniformRectPrimitiveUniformRooftops(rectPrimitive, quad);
rooftops.generate();

planner = rooftopMoM.core.greensMatrix.AssemblerChunkPlanner(elementCount);
context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);

assembler = rooftopMoM.core.greensMatrix.FullAssembler(planner);

outputConv = assembler.generateFull(context);

conjugateGradientSolver = rooftopMoM.core.solver.SpectrallyPreconditionedConjugateGradientGenome(context, outputConv);

tolerance = 1e-9;
maxIterCount = 2000;

outputConvSolver = conjugateGradientSolver.solve( ...
    'genome', genome, ...
    'excitation', excitation, ...
    'tolerance', tolerance, ...
    'maxIterCount', maxIterCount ...
);

fprintf('Number of the iterations: %d', outputConvSolver.count());
currentConv = outputConvSolver.postRunInfo().current

%
% Compare the currents.
%

normDifference = norm(currentNaive) - norm(currentConv)
