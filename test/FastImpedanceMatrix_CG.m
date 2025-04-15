clc; clearvars; format long;

%% Setup paths and parameters
currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 1;
countY = 1000;
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

%% Define phase factors
[M, N] = size(GT0xP);
[u, v] = ndgrid(0:(M - 1), 0:(N - 1));
Fx = exp((1i * pi * u) / gridSize(1));
FxStar = conj(Fx);
Fy = exp((1i * pi * v) / gridSize(2));
FyStar = conj(Fy);

%% Define excitation vector
feederPosition = ceil(countY / 2);
V = zeros(basisFunctions.getCount(), 1);
V(feederPosition) = 1;
B.x = [];  % Not used in this particular formulation.
B.y = V;

tolerance = 1e-6;
maxIterCount = 100000;

%% Run Conjugate Gradient solver using our FFT-based operator
tic;
[Jx, Jy, iter, res_norm, flag] = solve_CG(GT0xPF, GT0yPF, GP0PF, ...
    Fx, FxStar, Fy, FyStar, B, gridSize, tolerance, maxIterCount, M_pad, N_pad);
toc;

fprintf("CG norm = %d, iter count = %d\n", norm(Jy), iter);

%% --- Function Definitions ---

function [Jx, Jy, iter, res_norm, flag] = solve_CG(GT0xP, GT0yP, GP0P, ...
    Fx, FxStar, Fy, FyStar, B, gridSize, tol, maxIter, M_pad, N_pad)
    % Determine the sizes for Jx and Jy from the grid dimensions.
    nx = gridSize(1) - 1;  % Jx is of size (gridSize(1)-1) x gridSize(2)
    ny = gridSize(2);
    mx = gridSize(1);      % Jy is of size gridSize(1) x (gridSize(2)-1)
    my = gridSize(2) - 1;
    
    % Vectorize the excitation into a single column.
    % (B.x is empty in this setup; B.y provides the excitation).
    B_vec = [B.x(:); B.y(:)];
    
    % Define the operator G as a function handle.
    applyG = @(X) applyG_operator(X, GT0xP, GT0yP, GP0P, ...
        Fx, FxStar, Fy, FyStar, gridSize, M_pad, N_pad);
    
    % Use bicgstab to solve for J = [Jx; Jy]
    %[J_vec, flag, res_norm, iter] = bicgstab(applyG, B_vec, tol, maxIter);
    [J_vec, flag, res_norm, iter] = conjgrad(applyG, B_vec, tol, maxIter);

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
    
    % Reshape the input vector into its x and y components.
    Jx = reshape(J_vec(1:N_Jx), nx, ny);
    Jy = reshape(J_vec(N_Jx+1:end), mx, my);
    
    % Instead of using an expandJ function, use fft2’s ability to pad:
    JxPF = fft2(Jx, M_pad, N_pad);
    JyPF = fft2(Jy, M_pad, N_pad);
    
    % Combine FFT-domain multiplications to reduce number of FFTs.
    % For the x–component (Gx):
    Gxx_term = GT0xP + (1 - Fx) .* GP0P .* (1 - FxStar);
    GxxJx = ifft2( Gxx_term .* JxPF );
    GxyJx = ifft2( (1 - Fx) .* GP0P .* (1 - FyStar) .* JyPF );
    Gx = GxxJx + GxyJx;
    
    % For the y–component (Gy):
    GyxJy = ifft2( (1 - Fy) .* GP0P .* (1 - FxStar) .* JxPF );
    Gyy_term = GT0yP + (1 - Fy) .* GP0P .* (1 - FyStar);
    GyyJy = ifft2( Gyy_term .* JyPF );
    Gy = GyxJy + GyyJy;
    
    % Extract the portions corresponding to the original discretization.
    Gx_trim = Gx(1:nx, 1:ny);
    Gy_trim = Gy(1:mx, 1:my);
    
    % Concatenate the results into a single vector.
    GJ_vec = [Gx_trim(:); Gy_trim(:)];
end

function [x, flag, res_norm, iter] = conjgrad(applyA, b, tol, maxIter)
%CONJGRAD  Solve A*x = b by the Conjugate Gradient method,
%          where A is accessed via applyA(x) = A*x.
%
% Inputs:
%    applyA  - function handle that, given x, returns A*x
%    b       - right-hand side vector
%    tol     - tolerance for stopping criterion
%    maxIter - maximum number of iterations
%
% Outputs:
%    x        - approximate solution
%    flag     - 0 if converged, 1 otherwise
%    res_norm - final residual norm (relative to ||b||)
%    iter     - number of iterations taken

    % Initialize
    x = zeros(size(b));          % Start from x0 = 0 (can be changed if desired)
    r = b - applyA(x);           % r0 = b - A*x0
    p = r;                       % p0 = r0
    rTr = r.' * r;                % (r0, r0)
    bnorm = norm(b);             % Used to compute relative residual
    flag = 1;

    % CG iteration
    for iter = 1:maxIter
        Ap = applyA(p);          % A p_k
        pAp = p.' * Ap;           % (p_k, A p_k)
        
        % alpha_k = (r_k, r_k) / (p_k, A p_k)
        alpha = rTr / pAp;
        
        % x_{k+1} = x_k + alpha_k p_k
        x = x + alpha * p;
        
        % r_{k+1} = r_k - alpha_k A p_k
        r = r - alpha * Ap;
        
        % Check the new residual norm
        new_rTr = r.' * r;
        res_norm = sqrt(new_rTr) / bnorm;
        if res_norm < tol
            flag = 0;  % Converged
            break;
        end
        
        % beta_k = (r_{k+1}, r_{k+1}) / (r_k, r_k)
        beta = new_rTr / rTr;
        
        % p_{k+1} = r_{k+1} + beta_k p_k
        p = r + beta * p;
        
        % Prepare for next iteration
        rTr = new_rTr;
    end
end

