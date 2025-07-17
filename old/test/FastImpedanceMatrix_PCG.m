clc; clearvars; format long;

%% Setup paths and parameters
currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 100;
countY = 1;
sizeX = 1;
sizeY = 1;
gridSize = [countX, countY];
origin = [0, 0, 0];

order = 5;
lambda = 5;
k = 2 * pi / lambda;

%% Create mesh and related objects
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

%% Compute impedance matrix parts
tic;
[GT0x, GT0y] = strategyMag.compute(strategyCtx);
GP0 = strategyEl.compute(strategyCtx);
toc;

tic;
GT0xP = mom.FastImpedanceMatrix.CirculantMatrixEmbedding.compute(GT0x, gridSize);
GT0yP = mom.FastImpedanceMatrix.CirculantMatrixEmbedding.compute(GT0y, gridSize);
GP0P  = mom.FastImpedanceMatrix.CirculantMatrixEmbedding.compute(GP0, gridSize);
toc;

%% Define padded grid size for FFTs (automatic zero padding in fft2)
M_pad = 2 * gridSize(1);
N_pad = 2 * gridSize(2);

%% Precompute FFT of the circulant matrices
GT0xPF = fft2(GT0xP);
GT0yPF = fft2(GT0yP);
GP0PF  = fft2(GP0P);

%% Define phase factors (for the circulant embedding) 
[M, N] = size(GT0xP);
[u, v] = ndgrid(0:(M - 1), 0:(N - 1));
Fx = exp((1i * pi * u) / gridSize(1));
FxStar = conj(Fx);
Fy = exp((1i * pi * v) / gridSize(2));
FyStar = conj(Fy);

% For the (1,1) block, incorporate the phase factors symmetrically.
A11 = Fx .* GT0xPF .* FxStar + (1 - Fx) .* GP0PF .* (1 - FxStar);

% For the (2,2) block: use the corresponding phase factors.
A22 = Fy .* GT0yPF .* FyStar + (1 - Fy) .* GP0PF .* (1 - FyStar);

% For the off-diagonal blocks, the phase factors must be applied asymmetrically:
A12 = Fx .* GP0PF .* (1 - FyStar);
A21 = (1 - Fy) .* GP0PF .* FyStar;


%% Define excitation vector
feederPosition = ceil(countY / 2);
V = zeros(basisFunctions.getCount(), 1);
V(feederPosition) = 1;
B.x = [];  % Not used in this formulation.
B.y = V;

tolerance = 1e-6;
maxIterCount = 1000;

%% Run the preconditioned Conjugate Gradient solver
tic;
[Jx, Jy, iter, res_norm, flag] = solve_precond_CG(GT0xP, GT0yP, GP0P, ...
    Fx, FxStar, Fy, FyStar, A11, A12, A21, A22, B, gridSize, tolerance, maxIterCount, M_pad, N_pad);
toc;

fprintf("Preconditioned CG norm = %d, iter count = %d\n", norm(Jy), iter);

%% --- Function Definitions ---

function [Jx, Jy, iter, res_norm, flag] = solve_precond_CG(GT0xP, GT0yP, GP0P, ...
    Fx, FxStar, Fy, FyStar, A11, A12, A21, A22, B, gridSize, tol, maxIter, M_pad, N_pad)
    % Determine sizes for Jx and Jy
    nx = gridSize(1) - 1;  % Jx is (gridSize(1)-1) x gridSize(2)
    ny = gridSize(2);
    mx = gridSize(1);      % Jy is gridSize(1) x (gridSize(2)-1)
    my = gridSize(2) - 1;
    
    % Vectorize the excitation into a single column.
    % (B.x is empty; B.y provides the excitation)
    B_vec = [B.x(:); B.y(:)];
    
    % Define the forward operator as a function handle.
    applyG = @(X) applyG_operator(X, GT0xP, GT0yP, GP0P, ...
        Fx, FxStar, Fy, FyStar, gridSize, M_pad, N_pad);
    
    % Define the spectral preconditioner function.
    applyP = @(r) applySpectralPrecond(r, A11, A12, A21, A22, gridSize, M_pad, N_pad);
    
    % Use our preconditioned CG solver to solve for J = [Jx; Jy]
    [J_vec, flag, res_norm, iter] = precond_conjgrad(applyG, applyP, B_vec, tol, maxIter);
    
    % Extract Jx and Jy from the vectorized solution.
    Jx = reshape(J_vec(1:nx*ny), nx, ny);
    Jy = reshape(J_vec(nx*ny+1:end), mx, my);
end

function GJ_vec = applyG_operator(J_vec, GT0xP, GT0yP, GP0P, ...
        Fx, FxStar, Fy, FyStar, gridSize, M_pad, N_pad)
    % Recover sizes for splitting the vector.
    nx = gridSize(1) - 1; ny = gridSize(2);
    mx = gridSize(1);   my = gridSize(2) - 1;
    N_Jx = nx * ny;
    
    % Reshape input vector into x and y components.
    Jx = reshape(J_vec(1:N_Jx), nx, ny);
    Jy = reshape(J_vec(N_Jx+1:end), mx, my);
    
    % Zero-pad using fft2's padding feature.
    JxPF = fft2(Jx, M_pad, N_pad);
    JyPF = fft2(Jy, M_pad, N_pad);
    
    % Compute Fourier-domain multiplications.
    % For the x–component:
    Gxx_term = GT0xP + (1 - Fx) .* GP0P .* (1 - FxStar);
    GxxJx = ifft2(Gxx_term .* JxPF);
    GxyJx = ifft2((1 - Fx) .* GP0P .* (1 - FyStar) .* JyPF);
    Gx = GxxJx + GxyJx;
    
    % For the y–component:
    GyxJy = ifft2((1 - Fy) .* GP0P .* (1 - FxStar) .* JxPF);
    Gyy_term = GT0yP + (1 - Fy) .* GP0P .* (1 - FyStar);
    GyyJy = ifft2(Gyy_term .* JyPF);
    Gy = GyxJy + GyyJy;
    
    % Extract portions corresponding to the original discretization.
    Gx_trim = Gx(1:nx, 1:ny);
    Gy_trim = Gy(1:mx, 1:my);
    
    % Concatenate into a single vector.
    GJ_vec = [Gx_trim(:); Gy_trim(:)];
end

function z = applySpectralPrecond(r, A11, A12, A21, A22, gridSize, M_pad, N_pad)
    % This function applies the spectral preconditioner.
    % It splits the residual into x and y parts, zero-pads to the FFT grid,
    % computes the FFT, solves the 2x2 system at each Fourier frequency,
    % inverts the FFT and trims back to the original vector sizes.
    
    % Determine sizes.
    nx = gridSize(1) - 1; ny = gridSize(2);  % for x component
    mx = gridSize(1);    my = gridSize(2) - 1;  % for y component
    
    % Split the residual vector.
    r1 = reshape(r(1:nx*ny), nx, ny);
    r2 = reshape(r(nx*ny+1:end), mx, my);
    
    % Zero-pad the residuals.
    R1F = fft2(r1, M_pad, N_pad);
    R2F = fft2(r2, M_pad, N_pad);
    
    % Solve the 2x2 system in the Fourier domain at every frequency.
    % For each (u,v) we need to compute:
    %   det = A11 .* A22 - A12 .* A21
    %   Z1 = (A22 .* R1F - A12 .* R2F) ./ det
    %   Z2 = (-A21 .* R1F + A11 .* R2F) ./ det
    det = A11 .* A22 - A12 .* A21;
    
    % To avoid division by zero, one could add a small epsilon to det.
    epsilon = 1e-14;
    det = det + epsilon;
    
    Z1 = (A22 .* R1F - A12 .* R2F) ./ det;
    Z2 = (-A21 .* R1F + A11 .* R2F) ./ det;
    
    % Inverse FFT to return to spatial domain.
    z1_full = ifft2(Z1);
    z2_full = ifft2(Z2);
    
    % Trim to the original sizes.
    z1 = z1_full(1:nx, 1:ny);
    z2 = z2_full(1:mx, 1:my);
    
    % Concatenate into a single vector.
    z = [z1(:); z2(:)];
end

function [x, flag, res_norm, iter] = precond_conjgrad(applyA, applyP, b, tol, maxIter)
% PRECOND_CONJGRAD  Solve A*x = b using the Preconditioned Conjugate Gradient method,
% where A is accessed via applyA and the preconditioner M^{-1} is applied using applyP.
%
% Inputs:
%    applyA  - function handle that returns A*x
%    applyP  - function handle that returns the preconditioned vector M^{-1}r
%    b       - right-hand side vector
%    tol     - tolerance for stopping criterion (relative residual)
%    maxIter - maximum number of iterations
%
% Outputs:
%    x       - approximate solution
%    flag    - 0 if converged, 1 otherwise
%    res_norm- final residual norm (relative to norm(b))
%    iter    - number of iterations taken

    % Initialization.
    x = zeros(size(b));  % Initial guess x0 = 0.
    r = b - applyA(x);   % Compute initial residual.
    z = applyP(r);       % Preconditioned residual.
    p = z;               % Initial search direction.
    rz_old = r.' * z;    % Inner product.
    b_norm = norm(b);
    
    flag = 1;
    for iter = 1:maxIter
        Ap = applyA(p);
        alpha = rz_old / (p.' * Ap);
        x = x + alpha * p;
        r = r - alpha * Ap;
        res_norm = norm(r) / b_norm;
        if res_norm < tol
            flag = 0;
            break;
        end
        z = applyP(r);
        rz_new = r.' * z;
        beta = rz_new / rz_old;
        p = z + beta * p;
        rz_old = rz_new;
    end
end

