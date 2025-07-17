clc; clearvars; format long;

%% Setup paths and parameters (adjust as needed)
currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

% Grid parameters (for demonstration, pick numbers that yield non-empty arrays)
% For example, let gridSize = [Mx, Ny] with Mx>1 and Ny>1
countX = 100;  
countY = 100;  
sizeX  = 1;
sizeY  = 1;
gridSize = [countX, countY];
origin = [0, 0, 0];

order = 5;
lambda = 5;
k = 2 * pi / lambda;

%% Create mesh and related objects (using your framework)
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
strategyEl  = mom.FastImpedanceMatrix.ComputeElPart();

%% Compute the (unembedded) impedance parts in the spatial domain
tic;
[GT0x, GT0y] = strategyMag.compute(strategyCtx);
GP0 = strategyEl.compute(strategyCtx);
toc;

%% Perform circulant embedding of the impedance kernels
tic;
GT0xP = mom.FastImpedanceMatrix.CirculantMatrixEmbedding.compute(GT0x, gridSize);
GT0yP = mom.FastImpedanceMatrix.CirculantMatrixEmbedding.compute(GT0y, gridSize);
GP0P  = mom.FastImpedanceMatrix.CirculantMatrixEmbedding.compute(GP0, gridSize);
toc;

%% Precompute FFT of the embedded kernels
% These represent the frequency response of the circulant-embedded kernels.
GT0xPF = fft2(GT0xP);
GT0yPF = fft2(GT0yP);
GP0PF  = fft2(GP0P);

%% Define FFT phase factors based on the grid
[M, N] = size(GT0xP);
[u, v] = ndgrid(0:(M - 1), 0:(N - 1));
Fx = exp((1i * pi * u) / gridSize(1));
FxStar = conj(Fx);
Fy = exp((1i * pi * v) / gridSize(2));
FyStar = conj(Fy);

%% Define padded sizes for FFTs
% We use a padding such that fft2 automatically zero-pads the input.
M_pad = 2 * gridSize(1);
N_pad = 2 * gridSize(2);

%% Determine the size of unknowns and preallocate impedance matrix
% The currents are split into two parts:
%   - Jx has size: (gridSize(1)-1) x gridSize(2)
%   - Jy has size: gridSize(1) x (gridSize(2)-1)
nx = gridSize(1) - 1; ny = gridSize(2);
mx = gridSize(1);   my = gridSize(2) - 1;

N_Jx = nx * ny;
N_Jy = mx * my;
N_total = N_Jx + N_Jy;

% Preallocate the impedance matrix G (each column is the response to a unit impulse)
G = zeros(N_total, N_total);

%% Create the impedance matrix by applying the operator to each canonical basis vector
%fprintf('Building impedance matrix column-by-column...\n');
tic;
for col = 1:N_total
    % Create a unit impulse for the current column
    e = zeros(N_total, 1);
    e(col) = 1;
    
    % Compute the corresponding column using the FFT-based operator
    G(:, col) = applyG_operator(e, GT0xPF, GT0yPF, GP0PF, ...
                    Fx, FxStar, Fy, FyStar, gridSize, M_pad, N_pad);
    
    if mod(col, round(N_total/10)) == 0  % simple progress indicator
        %fprintf('Finished column %d of %d\n', col, N_total);
    end
end
toc;

disp(size(G));

%fprintf('Impedance matrix built: size %d x %d\n', size(G,1), size(G,2));

%% (Optional) Save or further process G
% save('ImpedanceMatrix.mat', 'G');

%% --- Operator Function ---
function GJ_vec = applyG_operator(J_vec, GT0xPF, GT0yPF, GP0PF, ...
    Fx, FxStar, Fy, FyStar, gridSize, M_pad, N_pad)
    % This function applies the impedance operator to the input vector J_vec.
    % It assumes that the input current vector is partitioned as:
    %   J_vec = [Jx(:); Jy(:)]
    % where:
    %   - Jx has size (gridSize(1)-1) x gridSize(2)
    %   - Jy has size gridSize(1) x (gridSize(2)-1)
    
    % Determine sizes for each block:
    nx = gridSize(1) - 1; ny = gridSize(2);
    mx = gridSize(1);   my = gridSize(2) - 1;
    N_Jx = nx * ny;
    
    % Unpack the vector into its components:
    Jx = reshape(J_vec(1:N_Jx), nx, ny);
    Jy = reshape(J_vec(N_Jx+1:end), mx, my);
    
    % Compute FFTs of the currents with zero-padding:
    JxPF = fft2(Jx, M_pad, N_pad);
    JyPF = fft2(Jy, M_pad, N_pad);
    
    % Compute the x–component of the response:
    % The operator acting on Jx, plus the cross term from Jy.
    % (The following multiplications assume the impedance operator is
    % assembled via a circulant embedding of the underlying Green's function.)
    Gxx_term = GT0xPF + (1 - Fx) .* GP0PF .* (1 - FxStar);
    GxxJx = ifft2( Gxx_term .* JxPF );
    GxyJy = ifft2( (1 - Fx) .* GP0PF .* (1 - FyStar) .* JyPF );
    Gx = GxxJx + GxyJy;
    
    % Compute the y–component of the response:
    GyxJx = ifft2( (1 - Fy) .* GP0PF .* (1 - FxStar) .* JxPF );
    Gyy_term = GT0yPF + (1 - Fy) .* GP0PF .* (1 - FyStar);
    GyyJy = ifft2( Gyy_term .* JyPF );
    Gy = GyxJx + GyyJy;
    
    % Extract the components corresponding to the original (unpadded) grid:
    % Note that:
    %   - Gx is computed on an M_pad x N_pad grid; we only want its (1:nx, 1:ny) part.
    %   - Similarly for Gy, we extract (1:mx, 1:my).
    Gx_trim = Gx(1:nx, 1:ny);
    Gy_trim = Gy(1:mx, 1:my);
    
    % Pack the results back into a vector:
    GJ_vec = [Gx_trim(:); Gy_trim(:)];
end
