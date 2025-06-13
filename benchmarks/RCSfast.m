clc; clearvars; format long;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

%profile clear;
%profile -memory on;

% constants
c0 = 299792458; % m/s

%
% inputs
%

kaMin = 1e-2;

% this value will determine the radiator size
kaMax = 1000;

% number of samples
N = 60;

numberOfRooftopsPerOneWavelength = 10;
rectangleSideRatio = 2;

ka = logspace(log10(kaMin), log10(kaMax), N);

times = zeros(2, N);

default = struct( ...
    'RCS',  0,  ...
    'RCSNorm', 0, ...
    'lx', 0, ...
    'ly', 0, ...
    'area', 0 ...
);

values = repmat(default, N, 1);

% unknown threshold
unknownCountThreshold = 1e7;

% quadrature order
order = 3;

% memory per one chunk element-wise
memoryElementwise = 1e7;
memory_MB = memoryElementwise * 8 / 1e6;

% CG settings
tolerance = 1e-5;
maxIterationCount = 10000;
restart = 50;

% plane wave parameters
direction = [0, 0, -1];

% polarization parallel to the y axis
polarization = [0, 1, 0];

% RCS config
rHat = [0, 0, -1];

% intensity of the incidental wave
Ei = 1;

%
% preparation
%

kaP = 10;

lastCurrentVector = [];

for n = 1:N
    if ka(n) < kaP
        patchCountPerShorterSide = ceil((kaP * numberOfRooftopsPerOneWavelength) / (pi * sqrt(1 + floor(rectangleSideRatio^2))));
    else
        patchCountPerShorterSide = ceil((ka(n) * numberOfRooftopsPerOneWavelength) / (pi * sqrt(1 + floor(rectangleSideRatio^2))));
    end

    % the shorter side will be aligned with the x axis
    countX = patchCountPerShorterSide;
    countY = floor(rectangleSideRatio) * patchCountPerShorterSide;

    fprintf("Size of the radiator: %d x %d\n", countX, countY);

    gridSize = [countX, countY];

    % the patch dimensios stays 1 m x 1 m
    sizeX = 1;
    sizeY = 1;

    origin = [0, 0, 0];

    % shift to center the mesh (probably not neccessery)
    shift = [countX/2, countY/2, 0];

    % mesh generation
    rectPrimitive = mom.MeshRectPrimitive(countX, countY, sizeX, sizeY, []);
    rectMesh = mom.MeshRect(rectPrimitive, origin - shift, 0, 0, 0);
    meshCluster = mom.Mesh2DCluster();
    meshCluster.loadMesh(rectMesh.getAsMesh());

    % obtain the geometrical size of the radiator
    [lx, ly] = rectPrimitive.getMeshSize();
    area = double(lx * ly);

    values(n).lx = lx;
    values(n).ly = ly;
    values(n).area = area;

    quadrature = mom.Quadrature(order);
    basisFunctions = mom.BasisFunctions(meshCluster);

    unknownCount = basisFunctions.getCount();
    fprintf("Unknown count: %d\n", unknownCount);

    if unknownCount > unknownCountThreshold
        error("Unknown count threshold, reduce the kaMax.\n");
    end

    fastRooftops = mom.FastImpedanceMatrix.RooftopBasisFunctions(quadrature, gridSize);

    planeWaveExcitatiton = mom.ImpedanceMatrix.ComputePlaneWaveExcitation();

    fprintf("\n[[Computing the sample %d]]\n\n", n);

    % calculation of the lambda from the ka
    lambda = ((double(patchCountPerShorterSide) * pi) / ka(n)) * sqrt(1 + rectangleSideRatio^2);
    k = 2 * pi / lambda;
    fprintf("lambda = %f\n", lambda);

    fastStrategyCtx = mom.ImpedanceMatrix.StrategyContext(quadrature, basisFunctions, fastRooftops, k);

    conjugateGradientSolver = mom.FastImpedanceMatrix.ConjugateGradientSolver(fastStrategyCtx);

%
% operator generation
%

    tic;

    conjugateGradientSolver.generateOperators(memory_MB);

    timeToComputeOperators = toc;

    fprintf("Time to compute the Green's operator: %f\n", timeToComputeOperators);

    times(1, n) = timeToComputeOperators;

%
% excitation - plane wave perpendicular to the radiator
%

    V = planeWaveExcitatiton.compute(fastStrategyCtx, direction, polarization);
    V = V(:);

    excitationVectorX = V(1:basisFunctions.getCountInDirection('x'));
    excitationVectorY = V((basisFunctions.getCountInDirection('x') + 1):end);

%
% PCG-FFT solver
%

    fprintf("Starting the solver ...\n");

    tic;

    [currentX, currentY, converged, residueNorm, iterationCount] = conjugateGradientSolver.solveCurrentPCG(excitationVectorX, excitationVectorY, tolerance, maxIterationCount, lastCurrentVector);

    timeToSolveTheSystem = toc;

    clear excitationVectorX excitationVectorY;

    fprintf("Time to solve using PCG-FFT: %f\n", timeToSolveTheSystem);
    fprintf("Number of iterations: %s\n", mat2str(iterationCount));

    times(2, n) = timeToSolveTheSystem;

%
% RCS computation
%

    % reshape the current
    I = [currentX(:); currentY(:)];

    lastCurrentVector = I;

    rho = mom.Operators.ComputeRCS.computeFromExcitation(fastStrategyCtx, I, V, k, Ei);
    %rho = mom.Operators.ComputeRCS.computeInDirection(fastStrategyCtx, I, k, rHat, Ei);
    fprintf("RCS = %f\n", rho);
    values(n).RCS = rho;
    values(n).RCSNorm = rho / area^2;
end


%%

RCS_vec = [values.RCS];
RCSNorm_vec = [values.RCSNorm];
lx_vec = [values.lx];
area_vec = [values.area];

lambda_vec = (pi ./ ka) .* sqrt(1 + rectangleSideRatio^2) .* double(lx_vec) .* double(sizeX);

pastelRed   = [1.0, 0.5, 0.5];
pastelGreen = [0.5, 1.0, 0.5];
pastelBlue = [0.5, 0.5, 1.0];
pastelYellow = [0.5, 1.0, 1.0];

figure;
loglog(ka, RCS_vec, 'o-', 'LineWidth', 1.5, 'MarkerSize', 10, 'Color', pastelRed);
hold on;
loglog(ka, RCSNorm_vec, 'o-', 'LineWidth', 1.5, 'MarkerSize', 10, 'Color', pastelBlue);
hold on;
loglog(ka, 4 * pi * area_vec.^2 ./ lambda_vec.^2, 'x-', 'LineWidth', 1.5, 'MarkerSize', 10, 'Color', pastelGreen);
hold on;
loglog(ka, 4 * pi ./ lambda_vec.^2, 'x-', 'LineWidth', 1.5, 'MarkerSize', 10, 'Color', pastelYellow);
grid on;

% Labeling
xlabel('ka');
ylabel('RCS');

legend({'MoM', 'Analytical'}, 'Location', 'best');
