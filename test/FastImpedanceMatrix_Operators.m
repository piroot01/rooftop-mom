clc; clearvars; format long;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 1;
countY = 500;
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
%{
function [Jx, Jy, iter, res_norm, flag] = solve_CG(GT0xP, GT0yP, GP0P, ...
    Fx, FxStar, Fy, FyStar, B, gridSize, tol, maxIter)
    % solve_CG solves G*J = B using a preconditioned CGS iterative method.
    %
    % Preconditioning is applied via the left preconditioner M_inv,
    % approximated here as a Jacobi (diagonal) preconditioner.
    
    %--- Sizes for Jx and Jy:
    nx = gridSize(1) - 1; ny = gridSize(2);
    mx = gridSize(1);   my = gridSize(2) - 1;
    
    %--- Vectorize the RHS:
    B_vec = [B.x(:); B.y(:)];
    
    %--- Define the operator application function:
    applyG = @(X) applyG_operator(X, GT0xP, GT0yP, GP0P, ...
        Fx, FxStar, Fy, FyStar, gridSize);
    
    %--- Compute an approximate diagonal of the operator G.
    % As a simple approximation, we use the zero-frequency (DC) term in the FFT.
    %
    % For the circulant embedding, note that the Fourier mode corresponding to u=v=0
    % typically gives a good representation of the average diagonal value.
    %
    % Preconditioning for x–component:
    Mx_approx = GT0xP(1,1) + (1 - Fx(1,1)) * GP0P(1,1) * (1 - FxStar(1,1));
    % Preconditioning for y–component:
    My_approx = GT0yP(1,1) + (1 - Fy(1,1)) * GP0P(1,1) * (1 - FyStar(1,1));
    
    % Build a preconditioner vector for the full current [Jx; Jy]:
    % For each degree of freedom in Jx and Jy, we assign the same diagonal value.
    M_inv_diag = [ repmat(1/Mx_approx, nx*ny, 1); 
                   repmat(1/My_approx, mx*my, 1) ];
               
    % Define the left preconditioner operator:
    precond = @(x) M_inv_diag .* x;
    
    %--- Use MATLAB's preconditioned CGS method:
    % Here, we pass the preconditioner as the first left preconditioner argument.
    [J_vec, flag, res_norm, iter] = cgs(applyG, B_vec, tol, maxIter, precond);
    
    %--- Unvectorize the solution back to Jx and Jy:
    Jx = reshape(J_vec(1:nx*ny), nx, ny);
    Jy = reshape(J_vec(nx*ny+1:end), mx, my);
end
%}

function [Jx, Jy, iter, res_norm, flag] = solve_CG(GT0xP, GT0yP, GP0P, ...
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

    %--- Vectorize B similarly:
    B_vec = [B.x(:); B.y(:)];

    %--- Define the operator application as an anonymous function.
    % This helper function receives a vector "X" that represents [Jx; Jy]
    % and returns G*X (also as a vector).
    applyG = @(X) applyG_operator(X, GT0xP, GT0yP, GP0P, ...
        Fx, FxStar, Fy, FyStar, gridSize);

    [J_vec, flag, res_norm, iter] = cgs(applyG, B_vec, tol, maxIter);

    %--- Unvectorize the solution back to Jx and Jy:
    Jx = reshape(J_vec(1:nx*ny), nx, ny);
    Jy = reshape(J_vec(nx*ny+1:end), mx, my);
end

%{
%%%%%%%%%%%%%% GMRES

function [Jx, Jy, iter, res_norm, flag] = solve_GMRES(GT0xP, GT0yP, GP0P, Fx, FxStar, Fy, FyStar, B, gridSize, tol, maxIter)
    % Sizes
    nx = gridSize(1) - 1; ny = gridSize(2);
    mx = gridSize(1);   my = gridSize(2) - 1;
    B_vec = [B.x(:); B.y(:)];

    % Define operator and preconditioner
    applyG = @(X) applyG_operator(X, GT0xP, GT0yP, GP0P, Fx, FxStar, Fy, FyStar, gridSize);

    % Solve with restarted GMRES
    restart = 30; % Adjust based on memory
    [J_vec, flag, res_norm, iter] = gmres(applyG, B_vec, restart, tol, maxIter);

    % Extract Jx, Jy
    Jx = reshape(J_vec(1:nx*ny), nx, ny);
    Jy = reshape(J_vec(nx*ny+1:end), mx, my);
end

%}

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

feederPosition = ceil(countY / 2);
V = zeros(basisFunctions.getCount(), 1);
V(feederPosition) = 1;

B.x = [];
B.y = V;

tolerance = 1e-5;
maxIterCount = 100000;

tic;

[Jx, Jy, iter, res_norm, flag] = solve_CG(GT0xPF, GT0yPF, GP0PF, Fx, FxStar, Fy, FyStar, B, gridSize, tolerance, maxIterCount);

toc;

iter

norm(Jy)

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
