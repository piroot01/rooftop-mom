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

countX = 500;
countY = 500;
sizeX = 1;
sizeY = 1;

order = 3;
elementCount = 1e6;
k = pi / (100 * sizeY * countY);
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

greensMatrixApplyHandle = @(I) greensMatrixApplicator.apply(I);

% Generate the transformation matrix from the loop-tree basis into the rooftop.
[loopTreeTrans, loopTreeData] = rooftopMoM.core.solver.altBasis.generateLoopTreeToRooftopTransformation( ...
    countX, countY, sizeX, sizeY);

omega = physconst("lightSpeed") * k;

% Create a modified greens matrix applicator that operates on the loop-tree basis.
customApplyHandle = @(I) rooftopMoM.core.solver.altBasis.applyLoopTreeFreqNormBasisPerm( ...
    I, greensMatrixApplyHandle, loopTreeTrans, loopTreeData, omega);

% Instantiate the solver.
solver = rooftopMoM.core.solver.ConjugateGradient(customApplyHandle);
%solver = rooftopMoM.core.solver.ConjugateGradient(greensMatrixApplyHandle);

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
