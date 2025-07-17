clc; clearvars; format long;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

% only dipole antennas will make sense
countX = 50;
countY = 50;
gridSize = [countX, countY];
sizeX = 1;
sizeY = 1;
origin = [0, 0, 0];

order = 4;
lambda = 2 * countY * sizeY;
k = 2 * pi / lambda;

memoryElementwise = 1e5;
memory_MB = memoryElementwise * 8 / 1e6;

myRectMeshPrimitive1 = mom.MeshRectPrimitive(countX, countY, sizeX, sizeY, []);
myRectMesh1 = mom.MeshRect(myRectMeshPrimitive1, origin, 0, 0, 0);
myMeshCluster = mom.Mesh2DCluster();
myMeshCluster.loadMesh(myRectMesh1.getAsMesh());

quadrature = mom.Quadrature(order);
basisFunctions = mom.BasisFunctions(myMeshCluster);

rooftops = mom.RooftopBasisFunctions(myMeshCluster, basisFunctions.getFaces());
rooftops.generatePoints(quadrature);

fastRooftops = mom.FastImpedanceMatrix.RooftopBasisFunctions(quadrature, gridSize);
fastRooftops.generatePoints();

fastStrategyCtx = mom.ImpedanceMatrix.StrategyContext(quadrature, basisFunctions, fastRooftops, k);

conjugateGradientSolver = mom.FastImpedanceMatrix.ConjugateGradientSolver(fastStrategyCtx);

conjugateGradientSolver.generateOperators(memory_MB);

impedanceMatrix = conjugateGradientSolver.computeImpedanceMatrix1();
preconditionerMatrix = conjugateGradientSolver.computePreconditionerMatrix();

MZ = preconditionerMatrix * impedanceMatrix;

eigImpedanceMatrix = eig(impedanceMatrix);
eigPreconditioner = eig(MZ);

% condition numbers (2-norm)
cond_Z   = cond(impedanceMatrix);
cond_PZ  = cond(MZ);

% spectral radii
rho_Z    = max(abs(eigImpedanceMatrix));
rho_PZ   = max(abs(eigPreconditioner));

fprintf('\nCondition numbers:\n');
fprintf('cond(Z) = %e\n', cond_Z);
fprintf('cond(M^{-1} Z) = %e\n\n', cond_PZ);

fprintf('Spectral radii:\n');
fprintf('rho(Z) = %e\n', rho_Z);
fprintf('rho(M^{-1} Z) = %e\n\n', rho_PZ);

fprintf("Number of unknowns = %d\n\n", basisFunctions.getCount());

feederPosition = ceil(basisFunctions.getCount() / 2);
V = zeros(basisFunctions.getCount(), 1);
V(feederPosition) = 1;

%
% CG-FFT solver
%

excitationVectorX = [];
excitationVectorY = V;

restart = 100;
tolerance = 1e-6;
maxIterationCount = 10000;

% CG
tic;

[~, ~, ~, ~, iterationCount] = conjugateGradientSolver.solveCurrentCG(excitationVectorX, excitationVectorY, tolerance, maxIterationCount);

t1 = toc;

iterationBoundCG = ceil(1/2 * sqrt(cond_Z) * log(2 / tolerance));

fprintf("Time to solve using CG-FFT: %f\n", t1);
fprintf("Number of iterations: %s\n", mat2str(iterationCount));
fprintf("Iteration bound: %d\n\n", iterationBoundCG);

tic;

[~, ~, ~, ~, iterationCount] = conjugateGradientSolver.solveCurrentPCG(excitationVectorX, excitationVectorY, tolerance, maxIterationCount);

t2 = toc;

iterationBoundPCG = ceil(1/2 * sqrt(cond_PZ) * log(2 / tolerance));

fprintf("Time to solve using PCG-FFT: %f\n", t2);
fprintf("Number of iterations: %s\n", mat2str(iterationCount));
fprintf("Iteration bound: %d\n\n", iterationBoundPCG);

% figure;
% semilogx(abs(eigImpedanceMatrix), zeros(length(eigPreconditioner), 1), 'x', 'MarkerSize',8, 'LineWidth',1.5);
% hold on;
% semilogx(abs(eigPreconditioner), zeros(length(eigPreconditioner), 1), 'o', 'MarkerSize',8, 'LineWidth',1.5);
% hold on;
% plot([1 1], [-1 1], '-', 'LineWidth', 1, 'Color', 'red');

% fprintf("\n\n");
% for k = 1:length(eigPreconditioner)
%     fprintf("            (%f, 2)\n", abs(eigPreconditioner(k)));
% end
% 
% fprintf("\n\n");
% for k = 1:length(eigImpedanceMatrix)
%     fprintf("            (%f, 1)\n", abs(eigImpedanceMatrix(k)));
% end
