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

countX = 5;
countY = 5;
sizeX = 1;
sizeY = 1;

rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
rect = rooftopMoM.mesh.UniformRect(rectPrimitive);
mesh = rect.generateMesh();

order = 3;
elementCount = 1e7;
lambda = (sizeX * countX) / 2;
k = 2 * pi / lambda;
V = 1;

tolerance = 1e-6;
maxIterCount = 2000;

quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);

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

greensMatrixApplicator = rooftopMoM.core.greensMatrix.FFTAcceleratedApplicator(context, spectralGreensMatrix);
greensMatrixApplicator.generate();

%% Toeplitz preconditioner
preconditionerApplicator = rooftopMoM.core.greensMatrix.preconditioner.FFTAcceleratedToeplitz(context, spectralGreensMatrix);
preconditionerApplicator.generate();

conjugateGradientSolver = rooftopMoM.core.solver.PreconditionedConjugateGradient(greensMatrixApplicator, preconditionerApplicator);

outputConvSolver1 = conjugateGradientSolver.solve( ...
    'excitation', excitation, ...
    'tolerance', tolerance, ...
    'maxIterCount', maxIterCount ...
);

currentConv1 = outputConvSolver1.postRunInfo().current;

%% Strang preconditioner applicator
centerScale = 1.2;
preconditioner = rooftopMoM.core.greensMatrix.preconditioner.StrangMatrixData(context);
preconditioner.generate(greensMatrix, centerScale);

preconditionerApplicator = rooftopMoM.core.greensMatrix.preconditioner.FFTAcceleratedCirculant(context, preconditioner);

conjugateGradientSolver = rooftopMoM.core.solver.PreconditionedConjugateGradient(greensMatrixApplicator, preconditionerApplicator);

outputConvSolver2 = conjugateGradientSolver.solve( ...
    'excitation', excitation, ...
    'tolerance', tolerance, ...
    'maxIterCount', maxIterCount ...
);

currentConv2 = outputConvSolver2.postRunInfo().current;

%% Huckle preconditioner applicator
centerScale = 1.3;
preconditioner = rooftopMoM.core.greensMatrix.preconditioner.HuckleMatrixData(context);
preconditioner.generate(greensMatrix, 'full', centerScale);

preconditionerApplicator = rooftopMoM.core.greensMatrix.preconditioner.FFTAcceleratedCirculant(context, preconditioner);

conjugateGradientSolver = rooftopMoM.core.solver.PreconditionedConjugateGradient(greensMatrixApplicator, preconditionerApplicator);

outputConvSolver3 = conjugateGradientSolver.solve( ...
    'excitation', excitation, ...
    'tolerance', tolerance, ...
    'maxIterCount', maxIterCount ...
);

currentConv3 = outputConvSolver3.postRunInfo().current;

%% GMRES
matrixApply = @(x) preconditionerApplicator.apply(greensMatrixApplicator.apply(x));

opts.verbose = 2;

tic;
%[currentGMRES, resNormGMRES, timeGMRES] = external.fastgmres.fastgmres(matrixApply, preconditionerApplicator.apply(excitation), tolerance, maxIterCount, [], [], [], opts);
%[currentGMRES, ~, ~, ~, resNormGMRES] = gmres(matrixApply, preconditionerApplicator.apply(excitation), 100, tolerance, maxIterCount);
%[currentGMRES, ~, ~, ~, ~, resNormGMRES] = external.minresqlp.minresqlp(matrixApply, preconditionerApplicator.apply(excitation), tolerance, maxIterCount, [], [], [], [], 1);
tGMRES = toc;


%% Compare number of iterations
fprintf('Number of the iterations: %d, %d, %d\n', outputConvSolver1.count(), outputConvSolver2.count(), outputConvSolver3.count());
fprintf('Solving time: %g s, %g s, %g s, %g s, %g s\n', outputConvSolver1.postRunInfo().solvingTime, outputConvSolver2.postRunInfo().solvingTime, outputConvSolver3.postRunInfo().solvingTime, tGMRES);

%% Compare the currents.
norm(currentConv1 - currentConv2)
norm(currentConv1 - currentConv3)
%norm(currentConv1 - currentGMRES)

%% Plot residue norm convergence based on iterations
close all;

iterationInfoToep = outputConvSolver1.iterationInfo();
resNormToep = cellfun(@(c) c.residueNorm, iterationInfoToep, 'UniformOutput', true);

iterationInfoStrang = outputConvSolver2.iterationInfo();
timeStrang = cellfun(@(c) c.residueNorm, iterationInfoStrang, 'UniformOutput', true);

iterationInfoHuckle = outputConvSolver3.iterationInfo();
resNormHuckle = cellfun(@(c) c.residueNorm, iterationInfoHuckle, 'UniformOutput', true);

y1 = double(resNormToep(:));
y2 = double(timeStrang(:));
y3 = double(resNormHuckle(:));
%y4 = double(resNormGMRES(:));

x1 = 1:numel(y1);
x2 = 1:numel(y2);
x3 = 1:numel(y3);
%x4 = 1:numel(y4);

figure;
hold on;
semilogy(x1, y1, '-x', 'LineWidth', 2, 'MarkerSize', 8);
semilogy(x2, y2, '-o', 'LineWidth', 2, 'MarkerSize', 8);
semilogy(x3, y3, '-d', 'LineWidth', 2, 'MarkerSize', 8);
%semilogy(x4, y4, '-x', 'LineWidth', 2, 'MarkerSize', 8);
yline(tolerance, 'r--', 'LineWidth', 2);
hold off;

ax = gca;
ax.YScale = 'log';

grid on;
xlabel('Iteration');
ylabel('Residual norm');
legend('Toep','Strang','Huckle', 'GMRES', 'tolerace', 'Location','best');

%% Plot residue norm convergence based on solver time
timeToep = cellfun(@(c) c.solverTime, iterationInfoToep, 'UniformOutput', true);
timeStrang = cellfun(@(c) c.solverTime, iterationInfoStrang, 'UniformOutput', true);
timeHuckle = cellfun(@(c) c.solverTime, iterationInfoHuckle, 'UniformOutput', true);

t1 = double(timeToep(:));
t2 = double(timeStrang(:));
t3 = double(timeHuckle(:));
%t4 = double(timeGMRES(:));

figure;
hold on;
semilogy(t1, y1, '-x', 'LineWidth', 2, 'MarkerSize', 8);
semilogy(t2, y2, '-o', 'LineWidth', 2, 'MarkerSize', 8);
semilogy(t3, y3, '-d', 'LineWidth', 2, 'MarkerSize', 8);
%semilogy(t4, y4, '-x', 'LineWidth', 2, 'MarkerSize', 8);
yline(tolerance, 'r--', 'LineWidth', 2);
hold off;

ax = gca;
ax.YScale = 'log';

grid on;
xlabel('Time [s]');
ylabel('Residual norm');
legend('Toep','Strang','Huckle', 'GMRES', 'tolerace', 'Location','best');