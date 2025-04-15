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

    [J_vec, flag, res_norm, iter] = bicgstab(applyG, B_vec, tol, maxIter);

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

tolerance = 1e-6;
maxIterCount = 100000;

% CG

tic;

[Jx, Jy, iter, res_norm, flag] = solve_CG(GT0xPF, GT0yPF, GP0PF, Fx, FxStar, Fy, FyStar, B, gridSize, tolerance, maxIterCount);

toc;

fprintf("CG norm = %d, iter count = %d\n", norm(Jy), iter);

% BiCG

function [Jx, Jy, iter, res_norm, flag] = solve_BiCG(GT0xP, GT0yP, GP0P, ...
    Fx, FxStar, Fy, FyStar, B, gridSize, tol, maxIter)
    % Determine sizes for Jx and Jy:
    nx = gridSize(1) - 1; ny = gridSize(2);
    mx = gridSize(1);   my = gridSize(2) - 1;

    % Vectorize the right-hand side:
    B_vec = [B.x(:); B.y(:)];

    % Nested function to define the operator. Note that we now pass GT0yP correctly.
    function y = applyG(x, flagStr, GT0xP, GT0yP, GP0P, Fx, FxStar, Fy, FyStar, gridSize)
        if strcmp(flagStr, 'notransp')
            y = applyG_operator(x, GT0xP, GT0yP, GP0P, Fx, FxStar, Fy, FyStar, gridSize);
        elseif strcmp(flagStr, 'transp')
            % Pass GT0yP as the second parameter!
            y = applyGConjugated_operator(x, GT0xP, GT0yP, GP0P, Fx, FxStar, Fy, FyStar, gridSize);
        end
    end

    % Create the function handle for bicg:
    applyG_handle = @(x, flagStr) applyG(x, flagStr, GT0xP, GT0yP, GP0P, Fx, FxStar, Fy, FyStar, gridSize);

    % Call bicg with the provided tolerance and maximum iterations:
    [J_vec, flag, res_norm, iter] = bicg(applyG_handle, B_vec, tol, maxIter);

    % Reshape the solution vector back into Jx and Jy:
    Jx = reshape(J_vec(1:nx*ny), nx, ny);
    Jy = reshape(J_vec(nx*ny+1:end), mx, my);
end

%tic;

%[Jx, Jy, iter, res_norm, flag] = solve_BiCG(GT0xPF, GT0yPF, GP0PF, Fx, FxStar, Fy, FyStar, B, gridSize, tolerance, maxIterCount);

%toc;

%fprintf("BiCG norm = %d, iter count = %d\n", norm(Jy), iter);

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

function GJ_vec = applyGConjugated_operator(J_vec, GT0xP, GT0yP, GP0P, ...
    Fx, FxStar, Fy, FyStar, gridSize)
    % Sizes for Jx and Jy:
    nx = gridSize(1) - 1; 
    ny = gridSize(2);
    mx = gridSize(1);   
    my = gridSize(2) - 1;
    N_Jx = nx * ny;
    
    % Extract the current components:
    Jx = reshape(J_vec(1:N_Jx), nx, ny);
    Jy = reshape(J_vec(N_Jx+1:end), mx, my);
    
    % Expand (pad) the currents as done in the forward operator:
    JxP = expandJ(Jx, 'x', gridSize);
    JyP = expandJ(Jy, 'y', gridSize);
    
    % Compute FFTs of the padded currents:
    JxPF = fft2(JxP);
    JyPF = fft2(JyP);
    
    % For the x–component:
    % Compute the “magnetic” term using the conjugated kernel:
    magTermF_adj = ifft2(conj(GT0xP) .* JxPF);
    % Compute the “electric” term:
    elTermF_adj = ifft2(conj((1 - Fx) .* GP0P .* (1 - FxStar)) .* JxPF);
    GxxJx_adj = magTermF_adj + elTermF_adj;
    
    % Compute the cross contribution from Jy (x–component):
    GxyJy_adj = ifft2(conj((1 - Fx) .* GP0P .* (1 - FyStar)) .* JyPF);
    
    % For the y–component:
    % Cross contribution from Jx:
    GyxJx_adj = ifft2(conj((1 - Fy) .* GP0P .* (1 - FxStar)) .* JxPF);
    % “Magnetic” term for y–component:
    magTermF_y_adj = ifft2(conj(GT0yP) .* JyPF);
    % “Electric” term for y–component:
    elTermF_y_adj = ifft2(conj((1 - Fy) .* GP0P .* (1 - FyStar)) .* JyPF);
    GyyJy_adj = magTermF_y_adj + elTermF_y_adj;
    
    % Sum the contributions to get the full adjoint applied current:
    Gx_adj = GxxJx_adj + GxyJy_adj;
    Gy_adj = GyxJx_adj + GyyJy_adj;
    
    % Since the FFTs yield values on the full padded grid,
    % trim the result to the original grid size:
    Gx_trim = Gx_adj(1:nx, 1:ny);
    Gy_trim = Gy_adj(1:mx, 1:my);
    
    % Combine back into a vector (ordering as in your forward operator):
    GJ_vec = [Gx_trim(:); Gy_trim(:)];
end

% GMRES precond

function [Jx, Jy, iter, relres, flag] = solve_GMRES(GT0xP, GT0yP, GP0P, ...
    Fx, FxStar, Fy, FyStar, B, gridSize, tol, maxIter, restart)
    % solve_GMRES_precond solves G*J = B (where J = [Jx; Jy] is the rooftop current)
    % using GMRES with a simple Jacobi preconditioner.
    %
    % Inputs:
    %   GT0xP, GT0yP, GP0P : FFT–domain matrices (size [M,N])
    %   Fx, FxStar, Fy, FyStar : FFT phase factors defined on a grid [M,N]
    %   B : structure with fields:
    %         B.x  -- excitation for the x–component (size: gridSize - [1,0])
    %         B.y  -- excitation for the y–component (size: gridSize - [0,1])
    %   gridSize : [M, N] full grid size used for padding/embedding
    %   tol : tolerance for the GMRES relative residual (e.g., 1e-6)
    %   maxIter : maximum number of GMRES iterations
    %   restart : restart parameter for GMRES (typical values: 30-50)
    %
    % Outputs:
    %   Jx, Jy : computed current components (shaped appropriately)
    %   iter : number of iterations performed (GMRES inner iterations may be nested)
    %   relres : final relative residual norm
    %   flag : convergence flag from GMRES
    %
    % Determine sizes for Jx and Jy:
    nx = gridSize(1) - 1;  % Jx has size (nx x gridSize(2))
    ny = gridSize(2);
    mx = gridSize(1);      % Jy has size (gridSize(1) x (gridSize(2)-1))
    my = gridSize(2) - 1;

    % Vectorize the excitation:
    B_vec = [B.x(:); B.y(:)];

    % Define the operator using your FFT-based application:
    applyG = @(X) applyG_operator(X, GT0xP, GT0yP, GP0P, ...
                Fx, FxStar, Fy, FyStar, gridSize);

    % Solve using GMRES with restart.
    % Note: gmres returns [solution, flag, relres, iter, residual_vector].
    [J_vec, flag, relres, iter, ~] = gmres(applyG, B_vec, restart, tol, maxIter);

    % Unvectorize the solution:
    Jx = reshape(J_vec(1:nx*ny), nx, ny);
    Jy = reshape(J_vec(nx*ny+1:end), mx, my);
end

%restart = 100;

%tic;

%[Jx, Jy, iter, res_norm, flag] = solve_GMRES(GT0xPF, GT0yPF, GP0PF, Fx, FxStar, Fy, FyStar, B, gridSize, tolerance, maxIterCount, restart);

%toc;

%fprintf("GMRES norm = %d, iter count = %s\n", norm(Jy), mat2str(iter));

