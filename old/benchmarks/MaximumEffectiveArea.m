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
kaMax = 700;

% number of samples
N = 50;

numberOfRooftopsPerOneWavelength = 10;
rectangleSideRatio = 2;

patchCountPerShorterSide = ceil((kaMax * numberOfRooftopsPerOneWavelength) / (pi * sqrt(1 + floor(rectangleSideRatio^2))));

% the main input of the proram
ka = logspace(log10(kaMin), log10(kaMax), N);

% unknown threshold
unknownCountThreshold = 1e7;

% quadrature order
order = 3;

% memory per one chunk element-wise
memoryElementwise = 1e7;
memory_MB = memoryElementwise * 8 / 1e6;

% CG settings
tolerance = 1e-4;
maxIterationCount = 10000;

% plane wave parameters
direction = [0, 0, -1];

% polarization parallel to the y axis
polarization = [0, 1, 0];

%
% preparation
%

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
myRectMeshPrimitive1 = mom.MeshRectPrimitive(countX, countY, sizeX, sizeY, []);
myRectMesh1 = mom.MeshRect(myRectMeshPrimitive1, origin - shift, 0, 0, 0);
myMeshCluster = mom.Mesh2DCluster();
myMeshCluster.loadMesh(myRectMesh1.getAsMesh());

quadrature = mom.Quadrature(order);
basisFunctions = mom.BasisFunctions(myMeshCluster);

unknownCount = basisFunctions.getCount();
fprintf("Unknown count: %d\n", unknownCount);

if unknownCount > unknownCountThreshold
    error("Unknown count threshold, reduce the kaMax.\n");
end

fastRooftops = mom.FastImpedanceMatrix.RooftopBasisFunctions(quadrature, gridSize);

planeWaveExcitatiton = mom.ImpedanceMatrix.ComputePlaneWaveExcitation();

times = zeros(2, N);

default = struct( ...
    'Aeff',      0, ...
    'AeffNorm',  0  ...
);

values = repmat(default, N, 1);

for n = 1:N
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

    [currentX, currentY, converged, residueNorm, iterationCount] = conjugateGradientSolver.solveCurrentPCG(excitationVectorX, excitationVectorY, tolerance, maxIterationCount);

    timeToSolveTheSystem = toc;

    fprintf("Time to solve using PCG-FFT: %f\n", timeToSolveTheSystem);
    fprintf("Number of iterations: %s\n", mat2str(iterationCount));

    times(2, n) = timeToSolveTheSystem;

%
% Effective area computation
%

    % reshape the current
    I = [currentX(:); currentY(:)];

    v = V(:)' * I(:);
    Aeff = (fastStrategyCtx.constants.Z0/4) * abs(v).^2 / real(v);
    fprintf("Aeff = %f\n", Aeff);
    values(n).Aeff = Aeff;

    AeffNorm = Aeff / (countX * countY * sizeX * sizeY);
    fprintf("AeffNorm = %f\n", AeffNorm);
    values(n).AeffNorm = AeffNorm;
end


%%

AeffNorm_vec = [values.AeffNorm];

figure;
loglog(ka, AeffNorm_vec, 'o-', 'LineWidth', 2, 'MarkerSize', 6, 'Color', 'red');
grid on;

% Labeling
xlabel('ka');
ylabel('A_{eff} / (l^2 / 2)');
