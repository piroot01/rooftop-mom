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

countX = 10;
countY = 20;
sizeX = 1;
sizeY = 1;

rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
rect = rooftopMoM.mesh.UniformRect(rectPrimitive);
mesh = rect.generateMesh();

order = 3;
elementCount = 1e7;
k = 2 * pi / (sizeY * countY);
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

% Greens matrix applicator: Z * I = V
greensMatrixApplicator = rooftopMoM.core.greensMatrix.FFTAcceleratedApplicator(context, spectralGreensMatrix);
greensMatrixApplicator.generate();

% Preconditioner applicator
preconditionerApplicator = rooftopMoM.core.greensMatrix.preconditioner.FFTAcceleratedToeplitz(context, spectralGreensMatrix);
preconditionerApplicator.generate();

n = greensMatrixApplicator.inputSize();
k = 6;
opts.tol = 1e-6;
opts.maxit = 1000;
%{
[~, vals] = eigs(@(x) greensMatrixApplicator.apply(x), n, k, 'largestabs', opts);
diag(vals)
noPrecondRho = max(abs(diag(vals)))
%}
[~, vals1] = eigs(@(x) preconditionerApplicator.apply(greensMatrixApplicator.apply(x)), n, n);
vals1abs = abs(diag(vals1));
toepPrecondRho = max(vals1abs)
minEig = min(vals1abs);

toepPrecondCond = toepPrecondRho / minEig

%% Strang preconditioner
centerScale = 1;
strangPreconditioner = rooftopMoM.core.greensMatrix.preconditioner.StrangMatrixData(context);
strangPreconditioner.generate(greensMatrix, centerScale);

preconditionerApplicator = rooftopMoM.core.greensMatrix.preconditioner.FFTAcceleratedCirculant(context, strangPreconditioner);

[~, vals2] = eigs(@(x) preconditionerApplicator.apply(greensMatrixApplicator.apply(x)), n, n);
vals2abs = abs(diag(vals2));
strangPrecondRho = max(vals2abs)
minEig = min(vals2abs);

strangPrecondCond = strangPrecondRho / minEig

%% T. Chan's preconditioner
tChanPreconditioner = rooftopMoM.core.greensMatrix.preconditioner.TChanMatrixData(context);
tChanPreconditioner.generate(greensMatrix);

preconditionerApplicator = rooftopMoM.core.greensMatrix.preconditioner.FFTAcceleratedCirculant(context, tChanPreconditioner);

[~, vals3] = eigs(@(x) preconditionerApplicator.apply(greensMatrixApplicator.apply(x)), n, n);
vals3abs = abs(diag(vals3));
tChanPrecondRho = max(vals3abs)
minEig = min(vals3abs);

tChanPrecondCond = tChanPrecondRho / minEig

%% R. Chan's preconditioner
rChanPreconditioner = rooftopMoM.core.greensMatrix.preconditioner.RChanMatrixData(context);
rChanPreconditioner.generate(greensMatrix);

preconditionerApplicator = rooftopMoM.core.greensMatrix.preconditioner.FFTAcceleratedCirculant(context, rChanPreconditioner);

[~, vals4] = eigs(@(x) preconditionerApplicator.apply(greensMatrixApplicator.apply(x)), n, n);
vals4abs = abs(diag(vals4));
rChanPrecondRho = max(vals4abs)
minEig = min(vals4abs);

rChanPrecondCond = rChanPrecondRho / minEig

%% Huckle preconditioner - half
centerScale = -1;
huckleHalfPreconditioner = rooftopMoM.core.greensMatrix.preconditioner.HuckleMatrixData(context);
huckleHalfPreconditioner.generate(greensMatrix, 'half', centerScale);

preconditionerApplicator = rooftopMoM.core.greensMatrix.preconditioner.FFTAcceleratedCirculant(context, huckleHalfPreconditioner);

[~, vals5] = eigs(@(x) preconditionerApplicator.apply(greensMatrixApplicator.apply(x)), n, n);
vals5abs = abs(diag(vals5));
huckleHalfPrecondRho = max(vals5abs)
minEig = min(vals5abs);

huckleHalfPrecondCond = huckleHalfPrecondRho / minEig

%% Huckle preconditioner - sqrt
huckleSqrtPreconditioner = rooftopMoM.core.greensMatrix.preconditioner.HuckleMatrixData(context);
huckleSqrtPreconditioner.generate(greensMatrix, 'sqrt', centerScale);

preconditionerApplicator = rooftopMoM.core.greensMatrix.preconditioner.FFTAcceleratedCirculant(context, huckleSqrtPreconditioner);

[~, vals6] = eigs(@(x) preconditionerApplicator.apply(greensMatrixApplicator.apply(x)), n, n);
vals6abs = abs(diag(vals6));
huckleSqrtPrecondRho = max(vals6abs)
minEig = min(vals6abs);

huckleSqrtPrecondCond = huckleSqrtPrecondRho / minEig

%% New preconditoner
preconditionerApplicator = rooftopMoM.core.greensMatrix.preconditioner.FFTAcceleratedToeplitz(context, spectralGreensMatrix);
preconditionerApplicator.generate();

function y = applyPreconditioner(greensMatrixApplicator, preconditionerApplicator, x)
    u = preconditionerApplicator.apply(x);

    w = greensMatrixApplicator.apply(u);

    w = 2 * x - w;

    y = preconditionerApplicator.apply(w);
end


[~, vals8] = eigs(@(x) applyPreconditioner(greensMatrixApplicator, preconditionerApplicator, greensMatrixApplicator.apply(x)), n, n);
vals8abs = abs(diag(vals8));
newPrecondRho = max(vals8abs)
minEig = min(vals8abs);

newPrecondCond = newPrecondRho / minEig

%% Green's matrix applicator spectra
[~, vals7] = eigs(@(x) greensMatrixApplicator.apply(x), n, n);
vals7abs = abs(diag(vals7));
greensMatrixApplicatorPrecondRho = max(vals7abs)
minEig = min(vals7abs);

greensMatrixApplicatorPrecondCond = greensMatrixApplicatorPrecondRho / minEig

all(real(diag(vals7)) > 0)


%% Plot the data
close all;

figure;
hold on;
%semilogx(vals7abs, 0 * ones(size(vals7abs)), '+', 'MarkerSize', 10, 'LineWidth', 2);
semilogx(vals1abs, ones(size(vals1abs)), '+', 'MarkerSize', 10, 'LineWidth', 2);
semilogx(vals2abs, 2 * ones(size(vals2abs)), '+', 'MarkerSize', 10, 'LineWidth', 2);
semilogx(vals3abs, 3 * ones(size(vals3abs)), '+', 'MarkerSize', 10, 'LineWidth', 2);
semilogx(vals4abs, 4 * ones(size(vals4abs)), '+', 'MarkerSize', 10, 'LineWidth', 2);
semilogx(vals5abs, 5 * ones(size(vals5abs)), '+', 'MarkerSize', 10, 'LineWidth', 2);
semilogx(vals6abs, 6 * ones(size(vals6abs)), '+', 'MarkerSize', 10, 'LineWidth', 2);
semilogx(vals8abs, 7 * ones(size(vals8abs)), '+', 'MarkerSize', 10, 'LineWidth', 2);
hold off;
ylim([-1; 8]);
xlim([-100, max(vals7abs + 100)]);
legend('Spectral Toeplitz preconditioner', 'Strang preconditioner', 'T. Chan preconditioner', 'R. Chan preconditioner', 'Huckle preconditioner - half', 'Huckle preconditioner - sqrt');

figure;
hold on;
plot(real(diag(vals1)), imag(diag(vals1)), 'o', 'MarkerSize', 10, 'LineWidth', 2);
plot(real(diag(vals2)), imag(diag(vals2)), 'x', 'MarkerSize', 10, 'LineWidth', 2);
plot(real(diag(vals6)), imag(diag(vals6)), 'D', 'MarkerSize', 10, 'LineWidth', 2);
plot(real(diag(vals7)), imag(diag(vals7)), 'd', 'MarkerSize', 10, 'LineWidth', 2);
hold off;
grid on;
legend('Spectral Toeplitz preconditioner', 'Strang preconditioner', 'Huckle preconditioner', 'Greens matrix');