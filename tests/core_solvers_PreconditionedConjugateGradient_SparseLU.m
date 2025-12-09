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

countX = 10000;
countY = 1;
sizeX = 1;
sizeY = 1;

order = 3;
elementCount = 1e6;
k = pi / (sizeX * countX);
V = 1;

quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);

rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);

%% Full impedance matrix assembler
%{
rect = rooftopMoM.mesh.UniformRect(rectPrimitive);
mesh = rect.generateMesh();

rooftops = rooftopMoM.basisFunctions.MeshUniformRooftops(mesh, quad);
rooftops.generate();

planner = rooftopMoM.core.impedanceMatrix.AssemblerChunkPlanner(elementCount);
context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);

fullAssembler = rooftopMoM.core.impedanceMatrix.FullAssembler(planner);

output = fullAssembler.generateFull(context);

Z = output.value();
%}

%% Convolutional approach
rooftops = rooftopMoM.basisFunctions.UniformRectPrimitiveUniformRooftops(rectPrimitive, quad);

rooftops.generate();

planner = rooftopMoM.core.greensMatrix.AssemblerChunkPlanner(elementCount);
context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);

assembler = rooftopMoM.core.greensMatrix.FullAssembler(planner);

greensMatrix = assembler.generateFull(context);

spectralGreensMatrix = rooftopMoM.core.greensMatrix.SpectralMatrixData(context);
spectralGreensMatrix.generate(greensMatrix);

greensMatrixApplicator = rooftopMoM.core.greensMatrix.FFTAcceleratedApplicator(context, spectralGreensMatrix);
greensMatrixApplicator.generate();

conjugateGradientSolver = rooftopMoM.core.solver.ConjugateGradient(greensMatrixApplicator);

mode = rooftopMoM.operator.utils.generateMiddlePortMode(rooftops.count());
excitation = rooftopMoM.operator.DeltaGapFeederExcitation.generate(rooftops.count(), mode, V);

tolerance = 1e-12;
maxIterCount = 5e2;

outputConvSolver1 = conjugateGradientSolver.solve( ...
    'excitation', excitation, ...
    'tolerance', tolerance, ...
    'maxIterCount', maxIterCount ...
);

fprintf('Number of the iterations: %d\n', outputConvSolver1.count());
currentConv = outputConvSolver1.postRunInfo().current;

% Preconditioned variant using sparse LU preconditioner
bandwidth = 3;
preconditioner = rooftopMoM.core.greensMatrix.preconditioner.SparseLU(context, bandwidth);

probeVector = preconditioner.generateProbeVector();
probeResponse = greensMatrixApplicator.apply(probeVector);
%probeResponse = Z * probeVector;
preconditioner.generate(probeResponse);

preconditionedSolver = rooftopMoM.core.solver.PreconditionedConjugateGradient(greensMatrixApplicator, preconditioner);

outputConvSolver2 = preconditionedSolver.solve( ...
    'excitation', excitation, ...
    'tolerance', tolerance, ...
    'maxIterCount', maxIterCount ...
);

fprintf('Number of the iterations: %d\n', outputConvSolver2.count());

%% GMRES

applyGreensMatrix = @(x) greensMatrixApplicator.apply(x);
applyPreconditioner = @(x) preconditioner.apply(x);

applyBoth = @(x) applyPreconditioner(applyGreensMatrix(x));

tic;
[currentGMRES, flag, relres, iter, resGMRES] = tfqmr(applyBoth, excitation, tolerance, maxIterCount);
toc;

%% Compare the convergences
close all;

iterationInfo = outputConvSolver1.iterationInfo();
resNorm = cellfun(@(c) c.residueNorm, iterationInfo, 'UniformOutput', true);

iterationInfoPrec = outputConvSolver2.iterationInfo();
resPrec = cellfun(@(c) c.residueNorm, iterationInfoPrec, 'UniformOutput', true);

y1 = double(resNorm(:));
y2 = double(resPrec(:));
y3 = double(resGMRES(:));

x1 = 1:numel(y1);
x2 = 1:numel(y2);
x3 = 1:numel(y3);

figure;
hold on;
semilogy(x1, y1, '-x', 'LineWidth', 2, 'MarkerSize', 8);
semilogy(x2, y2, '-o', 'LineWidth', 2, 'MarkerSize', 8);
semilogy(x3, y3, '-d', 'LineWidth', 2, 'MarkerSize', 8);
yline(tolerance, 'r--', 'LineWidth', 2);
hold off;

ax = gca;
ax.YScale = 'log';

grid on;
xlabel('Iteration');
ylabel('Residual norm');