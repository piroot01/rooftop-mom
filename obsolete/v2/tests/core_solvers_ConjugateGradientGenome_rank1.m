clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

loggerDir = fullfile(currentDir, '..', '3rd', 'matlab-logger', 'src');
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
countY = 1000;
sizeX = 1;
sizeY = 1;

rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
rect = rooftopMoM.mesh.UniformRect(rectPrimitive);
mesh = rect.generateMesh();

order = 5;
elementCount = 1e6;
k = pi / (sizeY * countY);
V = 1;

quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);

%
% solve the init current for a full antenna
%

rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
rooftops = rooftopMoM.basisFunctions.UniformRectPrimitiveUniformRooftops(rectPrimitive, quad);
rooftops.generate();

planner = rooftopMoM.core.greensMatrix.AssemblerChunkPlanner(elementCount);
context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);

assembler = rooftopMoM.core.greensMatrix.FullAssembler(planner);

outputConv = assembler.generateFull(context);

conjugateGradientSolver = rooftopMoM.core.solver.ConjugateGradientGenome(context, outputConv);

mode = rooftopMoM.operators.utils.generateMiddlePortMode(rooftops.count());
excitation = rooftopMoM.operators.DeltaGapFeederExcitation.generate(rooftops.count(), mode, V);

tolerance = 1e-9;
maxIterCount = 2000;
genome = true(rooftops.count(), 1);

outputConvSolver = conjugateGradientSolver.solve( ...
    'genome', genome, ...
    'excitation', excitation, ...
    'tolerance', tolerance, ...
    'maxIterCount', maxIterCount ...
);

fprintf('Number of the iterations: %d\n', outputConvSolver.count());
initCurrent = outputConvSolver.postRunInfo().current;

%
% rank1 perturbe the genome
%

%genome(randi([1, rooftops.count()])) = false;
%genome = genome | mode
genome = true(rooftops.count(), 1);
genome(end) = false;

outputConvSolver = conjugateGradientSolver.solve( ...
    'genome', genome, ...
    'excitation', excitation(genome), ...
    'tolerance', tolerance, ...
    'maxIterCount', maxIterCount ...
);

fprintf('Number of the iterations: %d\n', outputConvSolver.count());
currentPerturbed = outputConvSolver.postRunInfo().current;

%
% rank-1 perturbe using the woodbury
%

Irank1Perturbed = conjugateGradientSolver.rank1PerturbationRemoval(excitation, tolerance, maxIterCount, ~genome);
