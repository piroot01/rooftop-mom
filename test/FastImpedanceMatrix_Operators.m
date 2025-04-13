clc; clearvars; format long;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 1;
countY = 10;
sizeX = 1;
sizeY = 1;
gridSize = [countX, countY];
origin = [0, 0, 0];

order = 5;

lambda = 5;
k = 2 * pi / lambda;

myRectMeshPrimitive1 = mom.MeshRectPrimitive(countX, countY, sizeX, sizeY, []);
myRectMesh1 = mom.MeshRect(myRectMeshPrimitive1, origin, 0, 0, 0);
myMeshCluster = mom.Mesh2DCluster();
myMeshCluster.loadMesh(myRectMesh1.getAsMesh());

quadrature = mom.Quadrature(order);

basisFunctions = mom.BasisFunctions(myMeshCluster);

rooftopBasisFunctions = mom.FastImpedanceMatrix.RooftopBasisFunctions(quadrature, gridSize);
rooftopBasisFunctions.generatePoints();

strategyCtx = mom.ImpedanceMatrix.StrategyContext(quadrature, basisFunctions, rooftopBasisFunctions, k);
strategyMag = mom.FastImpedanceMatrix.ComputeMagPart();
strategyEl = mom.FastImpedanceMatrix.ComputeElPart();


tic;

[GT0x, GT0y] = strategyMag.compute(strategyCtx);
GP0 = strategyEl.compute(strategyCtx);

toc;

tic;

GT0xP = mom.FastImpedanceMatrix.CirculantMatrixEmbedding.compute(GT0x, gridSize);
GT0yP = mom.FastImpedanceMatrix.CirculantMatrixEmbedding.compute(GT0y, gridSize);

GP0P = mom.FastImpedanceMatrix.CirculantMatrixEmbedding.compute(GP0, gridSize);

toc;

JxSize = gridSize - [1, 0];
JySize = gridSize - [0, 1];

Jx = randn(JxSize);
Jy = randn(JySize);

function result = expandJ(J, direction, gridSize)
    if direction == 'x'
        M = gridSize(1) - 1;
        N = gridSize(2);

        result = zeros(2 * gridSize);
        result(1:M, 1:N) = J;
    else
        M = gridSize(1);
        N = gridSize(2) - 1;

        result = zeros(2 * gridSize);
        result(1:M, 1:N) = J;
    end
end

%{
JxP = expandJ(Jx, 'x', gridSize);
JyP = expandJ(Jy, 'y', gridSize);

[M, N] = size(GT0xP);
[u, v] = ndgrid(0:(M - 1), 0:(N - 1));

Fx = exp((1i * pi * u) / gridSize(1));
FxStar = conj(Fx);

Fy = exp((1i * pi * v) / gridSize(2));
FyStar = conj(Fy);

GT0xPF = fft2(GT0xP);
GT0yPF = fft2(GT0yP);
GP0PF = fft2(GP0P);
JxPF = fft2(JxP);
JyPF = fft2(JyP);

% compute GxxJx
magTermF = GT0xPF .* JxPF;
elTermF = (1 - Fx) .* GP0PF .* (1 - FxStar) .* JxPF;

magTerm = ifft2(magTermF);
elTerm = ifft2(elTermF);

GxxJx = magTerm + elTerm;

% compute GxyJy
elTermF = (1 - Fx) .* GP0PF .* (1 - FyStar) .* JyPF;

GxyJy = ifft2(elTermF);

% compute GyxJx
elTermF = (1 - Fy) .* GP0PF .* (1 - FxStar) .* JxPF;

GyxJx = ifft2(elTermF);

% compute GyyJy

magTermF = GT0yPF .* JyPF;
elTermF = (1 - Fy) .* GP0PF .* (1 - FyStar) .* JyPF;

magTerm = ifft2(magTermF);
elTerm = ifft2(elTermF);

GyyJy = magTerm + elTerm;

%}

function [Jx, Jy, iter, res_norm] = solve_CG(GT0xP, GT0yP, GP0P, ...
    Fx, FxStar, Fy, FyStar, B, gridSize, tol, maxIter)
    % solve_CG solves G*J = B for the unknown surface current J,
    % where J = [Jx; Jy] is composed of the x– and y–components of the rooftop
    % current, using a Conjugate Gradient iterative method.
    %
    % Inputs:
    %   GT0xP, GT0yP, GP0P : Precomputed FFT-domain matrices (size [M,N])
    %   Fx, FxStar, Fy, FyStar : FFT phase factors defined on a grid [M,N]
    %   B : structure with fields:
    %         B.x  -- excitation for the x–component (size: gridSize - [1,0])
    %         B.y  -- excitation for the y–component (size: gridSize - [0,1])
    %   gridSize : [M, N] size for the full grid (used for padding/embedding)
    %   tol : tolerance for the relative residual norm (e.g. 1e-6)
    %   maxIter : maximum number of CG iterations
    %
    % Outputs:
    %   Jx, Jy : computed current components (with sizes as described above)
    %   iter : number of iterations performed
    %   res_norm : final relative residual norm

    %--- Sizes for Jx and Jy based on your discretization:
    nx = gridSize(1) - 1; % Jx size: (gridSize(1)-1) x gridSize(2)
    ny = gridSize(2);
    mx = gridSize(1);
    my = gridSize(2) - 1; % Jy size: gridSize(1) x (gridSize(2)-1)

    %--- Initial guess: choose random initial currents (or zeros)
    if nx == 0
        Jx0 = nan(nx, ny);
    else
        Jx0 = randn(nx, ny);
    end

    if my == 0
        Jy0 = nan(mx, my);
    else
        Jy0 = randn(mx, my);
    end

    Jx0
    Jy0

    %--- Vectorize the unknowns into one column vector:
    J_vec = [Jx0(:); Jy0(:)];

    %--- Vectorize B similarly:
    B_vec = [B.x(:); B.y(:)]

    %--- Define the operator application as an anonymous function.
    % This helper function receives a vector "X" that represents [Jx; Jy]
    % and returns G*X (also as a vector).
    applyG = @(X) applyG_operator(X, GT0xP, GT0yP, GP0P, ...
        Fx, FxStar, Fy, FyStar, gridSize);

    %--- Compute initial residual:
    r = B_vec - applyG(J_vec)
    p = r;
    res_norm = norm(r) / norm(B_vec)
    iter = 0;
    fprintf('Iteration %d: relative residual = %e\n', iter, res_norm);

    %--- Main CG loop:
    while (res_norm > tol) && (iter < maxIter)
        Ap = applyG(p);
        alpha = (r' * r) / (p' * Ap);
        J_vec = J_vec + alpha * p;
        r_new = r - alpha * Ap;
        beta = (r_new' * r_new) / (r' * r);
        p = r_new + beta * p;
        r = r_new;
        iter = iter + 1;
        res_norm = norm(r) / norm(B_vec);
        fprintf('Iteration %d: relative residual = %e\n', iter, res_norm);
    end

    %--- Unvectorize the solution back to Jx and Jy:
    Jx = reshape(J_vec(1:nx*ny), nx, ny);
    Jy = reshape(J_vec(nx*ny+1:end), mx, my);
end

GT0xPF = fft2(GT0xP);
GT0yPF = fft2(GT0yP);
GP0PF = fft2(GP0P);

[M, N] = size(GT0xP);
[u, v] = ndgrid(0:(M - 1), 0:(N - 1));

Fx = exp((1i * pi * u) / gridSize(1));
FxStar = conj(Fx);

Fy = exp((1i * pi * v) / gridSize(2));
FyStar = conj(Fy);

% excitation vector
B.x = [];
B.y = [0, 0, 0, 0, 1, 0, 0, 0, 0];

tolerance = 1e-3
maxIterCount = 3;

solve_CG(GT0xPF, GT0yPF, GP0PF, Fx, FxStar, Fy, FyStar, B, gridSize, tolerance, maxIterCount)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GJ_vec = applyG_operator(J_vec, GT0xP, GT0yP, GP0P, ...
    Fx, FxStar, Fy, FyStar, gridSize)
    % applyG_operator applies the operator G (via FFTs) to the current vector J_vec.
    % This function implements the matrix–vector product using your FFT–based
    % formulation.
    %
    % The unknown current J has two parts:
    %   - Jx (size: (gridSize(1)-1) x gridSize(2))
    %   - Jy (size: gridSize(1) x (gridSize(2)-1))
    %
    % The operator G is defined via:
    %   Gx = Gxx * Jx + Gxy * Jy
    %   Gy = Gyx * Jx + Gyy * Jy
    %
    % where the individual terms are computed by FFTs as in your code.

    %--- Sizes for Jx and Jy:
    nx = gridSize(1) - 1; ny = gridSize(2);
    mx = gridSize(1);   my = gridSize(2) - 1;
    N_Jx = nx * ny;

    %--- Extract Jx and Jy from the vector:
    Jx = reshape(J_vec(1:N_Jx), nx, ny);
    Jy = reshape(J_vec(N_Jx+1:end), mx, my);

    %--- Embed Jx and Jy using your expansion function (assumed to be on the path):
    JxP = expandJ(Jx, 'x', gridSize);  % returns a 2*gridSize array
    JyP = expandJ(Jy, 'y', gridSize);

    %--- Compute FFTs of the embedded currents:
    JxPF = fft2(JxP);
    JyPF = fft2(JyP);

    %--- Compute the operator terms using FFT-domain multiplications:
    %  (i) For the x–component:
    %   GxxJx = ifft2( GT0xP .* JxPF ) + ifft2( (1-Fx) .* GP0P .* (1-FxStar) .* JxPF )
    magTermF = GT0xP .* JxPF;
    elTermF = (1 - Fx) .* GP0P .* (1 - FxStar) .* JxPF;
    GxxJx = ifft2(magTermF) + ifft2(elTermF);

    %   GxyJy = ifft2( (1-Fx) .* GP0P .* (1-FyStar) .* JyPF )
    elTermF = (1 - Fx) .* GP0P .* (1 - FyStar) .* JyPF;
    GxyJy = ifft2(elTermF);

    %  (ii) For the y–component:
    %   GyxJx = ifft2( (1-Fy) .* GP0P .* (1-FxStar) .* JxPF )
    elTermF = (1 - Fy) .* GP0P .* (1 - FxStar) .* JxPF;
    GyxJx = ifft2(elTermF);

    %   GyyJy = ifft2( GT0yP .* JyPF ) + ifft2( (1-Fy) .* GP0P .* (1-FyStar) .* JyPF )
    magTermF = GT0yP .* JyPF;
    elTermF = (1 - Fy) .* GP0P .* (1 - FyStar) .* JyPF;
    GyyJy = ifft2(magTermF) + ifft2(elTermF);

    %--- Combine the terms:
    % The full x–component (Gx) is the sum GxxJx + GxyJy.
    % The full y–component (Gy) is the sum GyxJx + GyyJy.
    Gx = GxxJx + GxyJy;
    Gy = GyxJx + GyyJy;

    %--- Since the FFTs produce values on the padded grid, extract the relevant
    % portions corresponding to the original rooftop discretization:
    Gx_trim = Gx(1:nx, 1:ny);
    Gy_trim = Gy(1:mx, 1:my);

    %--- Vectorize the results and concatenate:
    GJ_vec = [Gx_trim(:); Gy_trim(:)];
end
