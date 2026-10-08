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

countX = 300;
countY = 600;
sizeX = 1;
sizeY = 1;

rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
rect = rooftopMoM.mesh.UniformRect(rectPrimitive);
mesh = rect.generateMesh();

order = 3;
elementCount = 1e7;
lambda = (sizeY * countY) / 2;
k = 2 * pi / lambda;
V = 1;

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

% Greens matrix applicator
greensMatrixApplicator = rooftopMoM.core.greensMatrix.FFTAcceleratedApplicator(context, spectralGreensMatrix);
greensMatrixApplicator.generate();

% Strang preconditioner
strangPreconditioner = rooftopMoM.core.greensMatrix.preconditioner.StrangMatrixData(context);

% objective function
function iter = objectiveFunction(context, greensMatrix, greensMatrixApplicator, strangPreconditioner, excitation, centerScale)
    strangPreconditioner.generate(greensMatrix, centerScale);
    preconditionerApplicator = rooftopMoM.core.greensMatrix.preconditioner.FFTAcceleratedCirculant(context, strangPreconditioner);

    %{
    % Determine the spectral raduis and condition number
    n = greensMatrixApplicator.inputSize();
    [~, vals] = eigs(@(x) preconditionerApplicator.apply(greensMatrixApplicator.apply(x)), n, n);
    valsAbs = abs(diag(vals));
    rho = max(valsAbs);
    minEig = min(valsAbs);
    cond = rho / minEig;
    %}

    % Solve the system
    conjugateGradientSolver = rooftopMoM.core.solver.PreconditionedConjugateGradient(greensMatrixApplicator, preconditionerApplicator);

    tolerance = 1e-6;
    maxIterCount = 1e3;
    outputConvSolver = conjugateGradientSolver.solve( ...
        'excitation', excitation, ...
        'tolerance', tolerance, ...
        'maxIterCount', maxIterCount ...
    );

    iter = outputConvSolver.count()
end

objectiveFunctionHandle = @(c) objectiveFunction(context, greensMatrix, greensMatrixApplicator, strangPreconditioner, excitation, c);

opts = struct( ...
    'method',    'golden', ...
    'init_step', 0.1, ...     % starting step (choose something reasonable)
    'step_red',  0.5, ...
    'tol',       1e-6, ...    % smallest step (stopping tolerance)
    'maxIter',   1e3, ...
    'verbose',   true ...
);

[xg, fg, hg] = minimize1D(objectiveFunctionHandle, 1, opts);

%% Plot optimiser history
close all;

x_hist = hg.x(:);
f_hist = hg.f(:);
iters  = 1:numel(x_hist);

% Plot
figure; 
yyaxis left
plot(iters, x_hist, '-o','LineWidth',1);
ylabel('x value')

yyaxis right
plot(iters, f_hist, '-s','LineWidth',1);
ylabel('f(x)')

xlabel('Iteration')
title('Greedy Optimizer History')
grid on

%% Test the searched minima
iterDefault = objectiveFunctionHandle(1);
iterOpt = objectiveFunctionHandle(xg);

xg

improvement = iterDefault / iterOpt

%% Compare the spectra
close all;

function [vals, rho, cond] = getSpectra(context, greensMatrix, greensMatrixApplicator, hucklePreconditioner, centerScale)
    hucklePreconditioner.generate(greensMatrix, 'full', centerScale);
    preconditionerApplicator = rooftopMoM.core.greensMatrix.preconditioner.FFTAcceleratedCirculant(context, hucklePreconditioner);

    % Determine the spectral raduis and condition number
    n = greensMatrixApplicator.inputSize();
    [~, vals] = eigs(@(x) preconditionerApplicator.apply(greensMatrixApplicator.apply(x)), n, n);
    valsAbs = abs(diag(vals));
    rho = max(valsAbs);
    minEig = min(valsAbs);
    cond = rho / minEig;
end

[valsDefault, rhoDefault, condDefault] = getSpectra(context, greensMatrix, greensMatrixApplicator, hucklePreconditioner, 1);
[valsOpt, rhoOpt, condOpt] = getSpectra(context, greensMatrix, greensMatrixApplicator, hucklePreconditioner, xg);

fprintf('default: rho = %g, cond = %g\n', rhoDefault, condDefault);
fprintf('optimised: rho = %g, cond = %g\n', rhoOpt, condOpt);

valsDefaultAbs = abs(diag(valsDefault));
valsOptAbs = abs(diag(valsOpt));

figure;
hold on;
semilogx(valsDefaultAbs, 0 * ones(size(valsDefaultAbs)), '+', 'MarkerSize', 10, 'LineWidth', 2);
semilogx(valsOptAbs, ones(size(valsOptAbs)), '+', 'MarkerSize', 10, 'LineWidth', 2);
hold off;
ylim([-1 3]);
xlim([min(valsDefaultAbs) - 10, max(valsDefaultAbs) + 10]);
legend('default', 'optimised');

%% centerScale vs iter
close all;

centerScale = 1:0.01:1.2;
count = numel(centerScale);
iters = zeros(count, 1);

for n = 1:count
    iters(n) = objectiveFunctionHandle(centerScale(n));
end

figure;
plot(centerScale, iters);