clc; clearvars;

%profile clear;
%profile -memory on;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', '..', 'src');
addpath(srcDir);

loggerDir = fullfile(currentDir, '..','..', '3rd', 'matlab-logger', 'src');
addpath(loggerDir);

% Configure the logger.
logger.Logger.configure( ...
    'timeFormat', 'ss.SSS', ...
    'minLevel', logger.LogLevel.TRACE, ...
    'sinks', { ...
        logger.sink.ConsoleSink(), ...
    } ...
);

countX = 200000;
countY = 1;
sizeX = 1;
sizeY = 1;

order = 3;
elementCount = 5e6;
k = pi / (sizeX * countX);
V = 1;

tolerance = 1e-6;
maxIterCount = 1e3;

quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);

%
% Naive approach
%

rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
%{
rect = rooftopMoM.mesh.UniformRect(rectPrimitive);
mesh = rect.generateMesh();
rooftops = rooftopMoM.basisFunctions.MeshUniformRooftops(mesh, quad);
rooftops.generate();

planner = rooftopMoM.core.impedanceMatrix.AssemblerChunkPlanner(elementCount);

context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);

elAssembler = rooftopMoM.core.impedanceMatrix.ElectricPartAssembler(planner);

out = elAssembler.generateFull(context);
Ze = out.value();
%}
% Prep data for the iterative solver

rooftops = rooftopMoM.basisFunctions.UniformRectPrimitiveUniformRooftops(rectPrimitive, quad);
rooftops.generate();
planner = rooftopMoM.core.greensMatrix.AssemblerChunkPlanner(elementCount);
context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);
assembler = rooftopMoM.core.greensMatrix.FullAssembler(planner);
greensMatrix = assembler.generateFull(context);

mode = rooftopMoM.operator.utils.generateMiddlePortMode(rooftops.count());
excitation = rooftopMoM.operator.DeltaGapFeederExcitation.generate(rooftops.count(), mode, V);

%% Alternative electric Green's matrix application approach
% Embed electric part into a circulant
[cX, cY] = size(greensMatrix.electric);
electricCirc = rooftopMoM.core.greensMatrix.utils.embedToeplitz4ToCirculant(greensMatrix.electric, cX, cY);

% NOTE: Use only the first column as the fft runs only in the X direction.
electricCirc = electricCirc(:, 1);

% Build the standard shift operators in spectra
[circSizeX, circSizeY] = size(electricCirc);
[forwardX, backwardX, ~, ~] = rooftopMoM.core.greensMatrix.utils.generateSpectralShiftOperators(circSizeX, circSizeY);

% build the spectral greens matrix and its applicator
function gxxf = buildSpectralGreensMatrix(electricCirc, forwardX, backwardX)
    gxxf = (1 - forwardX) .* (1 - backwardX) .* fft(electricCirc);
end

function output = applySpectralGreensMatrix(input, gxxf)
    inputLen = length(input);
    circLen = length(gxxf);

    inputSpectr = fft(input, circLen);
    output = ifft(gxxf .* inputSpectr);
    output = output(1:inputLen);
end

% Test the solver via build in CG solver
gxxf = buildSpectralGreensMatrix(electricCirc, forwardX, backwardX);


applySpectralGreensMatrixHandle = @(x) applySpectralGreensMatrix(x, gxxf);

%{
t0 = tic;
[Icg, flag0, ~ ,iter0] = cgs(applySpectralGreensMatrixHandle, excitation, tolerance, maxIterCount);
solvingTime = toc(t0);
fprintf('Unpreconditioned variant: flag0 = %d, iteration count = %d, solving time = %g s.\n', flag0, iter0, solvingTime);

% Solve for current using the \ operator
%Iref = Ze \ excitation;
%}
% Compare the iterative and direct solver
%normDiff = norm(Iref - Icg);
%fprintf('Norm difference is %g.\n', normDiff);

%% New custom preconditioner for the electric part of the impedance matrix based on T. Chan's circulant estimator and discrete sine transform
% The idea is to split the preconditioner into two separate preconditioners and then we use the rule (AB)^(-1) = B^(-1)A^(-1), where A is the preconditioner for the shift operator
% The shift operator in time domain is a circulant matrix with the first column of [2, -1, 0, ...., 0, -1] such matrix is nasty (bad condition number and the spectra contains zeros
% which prevents us from using the actual circulant in fft in hadamard division based preconditioner) and it is possible as we see later to reduce the size of the preconditioner to
% NxN, because we know the exact form of the complementary Toeplitz matrix, whose generating column is [2, -1, 0, ..., 0], the eigenvalues of such matrix are 
% lambda_k = 2 - 2 * cos(k * pi / (N + 1)) and therefore the preconditioner T = S L S^(T), where S is the matrix of discrete sine transformation and L is diagonal matrix with lambda_k
% on its diagonal. As for the electric preconditioner we use T. Chan's preconditioner minimising the Frobeniusnorm between the E (generating side of the Toeplitz matrix) and C
% (the set of all circulants)

% First we compute the spectra of the shift operator
function lambda = generateShiftEigs(cX)
    lambdaCount = cX;
    lambda = 2 - 2 * cos((1:lambdaCount).' * pi ./ (lambdaCount + 1));
end

function y = dst1(x)
    n = numel(x);
    % build odd extension and use FFT
    X = fft([0; x(:); 0; -flipud(x(:))]);
    % indices: MATLAB 1-based, take elements 2..n+1, multiply by -1/2
    y = -0.5 * imag(X(2:n+1));
end

% inverse of the (unnormalized) DST-I is the same DST-I up to sign/scale,
% but we will handle the orthonormalization below.
function y = applyShiftPrecond(b, lambda)
    % lambda : vector of eigenvalues (length n)
    n = numel(b);
    % normalization factor for orthonormal sine basis: sqrt(2/(n+1))
    alpha = sqrt(2/(n+1));
    % forward transform (S' * b). For our symmetry DST-I gives S'*b = S*b
    scoeff = alpha * dst1(b);        % = S_mat' * b
    % diagonal solve in spectral domain
    scoeff = scoeff ./ lambda(:);
    % inverse transform
    y = alpha * dst1(scoeff);        % = S_mat * (Lambda^{-1} * S_mat' * b)
end

% Applicator of the T. Chan's preconditioner
function output = applyTChanPrecond(input, lambda)
    inputLen = length(input);
    lambdaCount = numel(lambda);
    output = ifft(fft(input, lambdaCount) ./ lambda);
    output = output(1:inputLen);
end

function output = applyToeplitzPrecond(input, gxxf)
    inputLen = length(input);
    gxxfLen = length(gxxf);
    output = ifft(fft(input, gxxfLen) ./ gxxf);
    output = output(1:inputLen);
end

% Now lets build the both preconditioners and theirs handle as well
% First the shift operator
%lambdaShiftOperator = generateShiftEigs(cX);
%applyShiftPrecondHandle = @(x) applyShiftPrecond(x, lambdaShiftOperator);

% Then the T. Chan's preconditioner of the electric part
%tChan = rooftopMoM.core.greensMatrix.preconditioner.utils.generateTChanPreconditionerKernel(gxxf);
%lambdaTChan = fft(tChan);
%[forwardX, backwardX, ~, ~] = rooftopMoM.core.greensMatrix.utils.generateSpectralShiftOperators(length(tChan), 1)
%lambdaTChan = buildSpectralGreensMatrix(tChan, forwardX, backwardX)

%applyShiftPrecondHandle = @(x) applyShiftPrecond(x, lambda);

% Define the Combined Preconditioner
% This function represents M^(-1). It applies the shift inverse, then the electric inverse.
% Mathematically: P_inv = P_electric_inv * P_shift_inv
%applyToeplitzPrecondHandle = @(x) applyToeplitzPrecond(x, fft(electricCirc));
%applyBothPrecondHandle = @(x) applyToeplitzPrecondHandle(applyShiftPrecondHandle(x));
%{
applyTChanPrecondHandle = @(x) applyTChanPrecond(x, lambdaTChan);

% Solve using CGS with the Preconditioner argument
t0 = tic;
[Icg_prec1, flag1, relres, iter1] = cgs(applySpectralGreensMatrixHandle, excitation, tolerance, maxIterCount, applyTChanPrecondHandle);
solvingTime = toc(t0);
fprintf('Preconditioned variant: flag0 = %d, iteration count = %d, solving time = %g s.\n', flag1, iter1, solvingTime);

% Compare the iterative and direct solver
normDiff = norm(Icg - Icg_prec1);
fprintf('Norm difference is %g.\n', normDiff);
%}
% Generate the preconditioner from the ifft(gxx)

G = greensMatrix.electric; 
[cX, cY] = size(G);

% 2. Allocate output
% We perform the calculation directly in the spatial domain.
gxx_half = zeros(cX, cY, 'like', G); 

if cX > 1
    % --- Internal Points ---
    % Standard 2*x[n] - x[n-1] - x[n+1]
    gxx_half(2:end-1, :) = 2*G(2:end-1, :) - G(1:end-2, :) - G(3:end, :);

    % --- Boundary Points (Reflective/Neumann) ---
    % Top: Neighbor is index 2. -> 2*x[1] - x[2] - x[2] = 2*(x[1] - x[2])
    gxx_half(1, :) = 2 * (G(1, :) - G(2, :));

    % Bottom: Neighbor is index end-1. -> 2*x[end] - x[end-1] - x[end-1]
    gxx_half(end, :) = 2 * (G(end, :) - G(end-1, :));
else
    % Edge case: 1xN input. Laplacian of a single point depends on physics, 
    % but with this reflection logic, 2*x[1] - x[1] - x[1] = 0.
    gxx_half = zeros(cX, cY, 'like', G); 
end

tChan = rooftopMoM.core.greensMatrix.preconditioner.utils.generateTChanPreconditionerKernel(gxx_half);
lambdaTChan = fft(tChan);

applyTChanPrecondHandle = @(x) applyTChanPrecond(x, lambdaTChan);

% Solve using CGS with the Preconditioner argument
t0 = tic;
[Icg_gxx_half, flag1, relres, iter1] = cgs(applySpectralGreensMatrixHandle, excitation, tolerance, maxIterCount, applyTChanPrecondHandle);
solvingTime = toc(t0);
fprintf('Preconditioned variant gxx half: flag0 = %d, iteration count = %d, solving time = %g s.\n', flag1, iter1, solvingTime);

%% Generalized Black-Box "Probe" Preconditioner
% Automatically extracts a banded Toeplitz approximation of arbitrary width.

% USER PARAMETER: Control the width of the preconditioner here.
% bandwidth = 1  -> Tridiagonal (Size 3 stencil: captures neighbors)
% bandwidth = 5  -> Banded (Size 11 stencil: captures wider interactions)
bandwidth = 10; 

fprintf('Probing operator with bandwidth = %d...\n', bandwidth);

% 1. Create a Probe Vector (Dirac Delta in the middle)
N = length(excitation);
midNode = floor(N / 2);

% Safety check: Ensure mesh is large enough for the requested bandwidth
if (midNode - bandwidth < 1) || (midNode + bandwidth > N)
    error('Mesh size N=%d is too small for bandwidth=%d. Increase countX or reduce bandwidth.', N, bandwidth);
end

probeVector = zeros(N, 1);
probeVector(midNode) = 1;

% 2. Apply the Fast Operator (Costs 1 FFT)
% This recovers the central column/row of the implicit matrix.
probeResponse = applySpectralGreensMatrixHandle(probeVector);

% 3. Generalized Extraction Logic
% We extract the window of values centered around the probe node.
% Indices: e.g., for BW=1 -> [mid-1, mid, mid+1]
extractionIndices = (midNode - bandwidth) : (midNode + bandwidth);
extractedValues = probeResponse(extractionIndices);

% 4. Construct the Sparse Banded Matrix
% Map the extracted values to diagonals. 
% Note on Indexing:
% A value at index (mid - k) corresponds to the Super-Diagonal (+k).
% A value at index (mid + k) corresponds to the Sub-Diagonal (-k).
% Therefore, the diagonal indices corresponding to our extracted window
% are: [+bandwidth, ..., 0, ..., -bandwidth].

diagonals = bandwidth : -1 : -bandwidth;

% Create the B matrix for spdiags.
% spdiags expects B to have columns corresponding to the diagonals in 'd'.
% Since we assume a Toeplitz system (translation invariant), we replicate
% the extracted values down the columns.
B = repmat(extractedValues(:).', N, 1);

% Build the sparse matrix
P_probe = spdiags(B, diagonals, N, N);

% 5. Factorize Once (LU Decomposition)
% As bandwidth increases, this becomes slightly more expensive, 
% but it is still O(N) for small bandwidths.
[L_probe, U_probe] = lu(P_probe);
applyProbePrecond = @(x) U_probe \ (L_probe \ x);

%% LU preconditioner based on the direct assembly
[P, x0] = create_preconditioner(greensMatrix.magneticXX, greensMatrix.electric, 1);

probeVector = zeros(cX - 1, 1);
probeVector(1) = 1;
x2 = applySpectralGreensMatrixHandle(probeVector);

% ---------------------------------------------------------
% Solve using CGS with the Learned Preconditioner
% ---------------------------------------------------------
fprintf('Starting CGS with Bandwidth-%d Preconditioner...\n', bandwidth);
t0 = tic;

[Icg_prec2, flag, relres, iter] = cgs(...
    applySpectralGreensMatrixHandle, ... 
    excitation, ...
    tolerance, ...
    maxIterCount, ...
    applyProbePrecond ... 
);

solveTime = toc(t0);

fprintf('--------------------------------------------------\n');
fprintf('Probe Results (Bandwidth %d):\n', bandwidth);
fprintf('Flag: %d (0=converged)\n', flag);
fprintf('Iterations: %d\n', iter);
fprintf('Time: %g s\n', solveTime);
fprintf('Relative Residual: %g\n', relres);
fprintf('--------------------------------------------------\n');

%{
%% Custom PCG solver
% 1. Define your Combined Operator (S * E)
A_fun = @(x) applySpectralGreensMatrixHandle(x); 

% 2. Define your Chained Preconditioner (M2_inv * M1_inv)
% This ensures M_fun applies the exact logic: M^-1 * r
M_fun = @(r) applyToeplitzPrecondHandle(applyShiftPrecondHandle(r));

% 3. Run the solver
[I_sol, info] = solve_pcg_standalone(A_fun, M_fun, zeros(size(excitation)), excitation, 1e-6, 1000);

% 4. Check results
if info.converged
    disp(['Converged in ' num2str(info.iterations) ' iterations.']);
else
    disp('Did not converge.');
end
%}
function [current, output] = solve_pcg_standalone(A_fun, M_fun, initCurrent, excitation, tolerance, maxIterCount)
    % solve_pcg_standalone: Core Spectrally Preconditioned Conjugate-Gradient loop.
    %
    %   Inputs:
    %       A_fun        : Function handle for Matrix operator (A * x)
    %       M_fun        : Function handle for Preconditioner (M^-1 * x)
    %       initCurrent  : Initial guess column vector.
    %       excitation   : Right-hand side column vector.
    %       tolerance    : Relative tolerance.
    %       maxIterCount : Maximum iterations.
    %
    %   Returns:
    %       current      : The final solution vector.
    %       output       : Struct containing convergence history and stats.

    % Initialize output struct
    output = struct();
    output.iterationInfo = struct('residueNorm', [], 'solverTime', []);
    
    current = initCurrent;
    
    % Calculate initial residual: r = b - A*x
    residue = excitation - A_fun(current);
    
    % Apply preconditioner: z = M^-1 * r
    preconditionedResidue = M_fun(residue);
    
    direction = preconditionedResidue;
    
    % Note: Using .' (transpose) instead of ' (conjugate transpose) 
    % preserves support for Complex Symmetric systems (COCG) common in MoM.
    rTr = residue.' * preconditionedResidue;
    
    excitationNorm = norm(excitation);
    converged = false;
    
    t0 = tic;

    for iteration = 1:maxIterCount
        t1 = tic;

        gDirection = A_fun(direction);
        
        alpha = rTr / (direction.' * gDirection);
        
        current = current + alpha * direction;
        residue = residue - alpha * gDirection;
        
        preconditionedResidue = M_fun(residue);
        
        rTrNew = residue.' * preconditionedResidue;
        residueNorm = abs(sqrt(rTrNew) / excitationNorm);

        beta = rTrNew / rTr;
        direction = preconditionedResidue + beta * direction;
        rTr = rTrNew;

        % Insert iteration info
        iterationTime = toc(t1);
        
        % Store history in vectors for easy plotting later
        output.iterationInfo.residueNorm(end+1) = residueNorm;
        output.iterationInfo.solverTime(end+1) = toc(t0);

        if residueNorm < tolerance
            converged = true;
            break;
        end
    end

    solvingTime = toc(t0);

    % Populate final post-run info
    output.converged = converged;
    output.solvingTime = solvingTime;
    output.finalResidueNorm = residueNorm;
    output.iterations = iteration;
end

function [P, h] = create_preconditioner(M, E, lx)
    % Inputs:
    % M: Column vector of size (N-1) x 1
    % E: Column vector of size N x 1
    % lx: Scalar value
    
    % 1. Determine dimensions
    % The system size is determined by M (N-1)
    K = length(M); % This is N-1 in the text's notation
    
    % 2. Prepare the component vectors
    % Calculate E_diff term: E(m) - E(m+1)
    % Note: E is size K+1, so we take diff of first K vs second K elements
    E_diff = E(1:end-1) - E(2:end);
    
    % Construct the difference kernel 'd' for x(m) - x(m-1)
    % This corresponds to [1, -1, 0, ...] in the time domain
    d = zeros(K, 1);
    d(1) = 1;
    if K > 1
        d(2) = -1;
    end
    
    % 3. Compute the Impulse Response 'h' (First column of Z)
    % We use FFT to perform the circulant convolution efficiently.
    % Formula: h = M + (1/lx^2) * (E_diff circ_conv d)
    
    % Compute convolution in spectral domain
    conv_term_f = fft(E_diff) .* fft(d);
    conv_term = ifft(conv_term_f);
    
    % Combine to get the full first column of Z
    h = M + (1/lx^2) * conv_term;
    
    % 4. Extract Diagonals for the Preconditioner
    % Z is circulant, so diagonals are constant.
    % Z_ij = h(i-j). 
    % Main diagonal (i=j): h(0) -> Matlab index 1
    % Sub diagonal (i=j+1): h(1) -> Matlab index 2
    % Super diagonal (i=j-1): h(-1) -> Wraps to Matlab index end
    
    val_main = h(1);
    val_sub  = h(2);
    val_sup  = h(end);
    
    % 5. Construct Sparse Band Matrix
    % We create a tridiagonal matrix.
    % spdiags takes columns of values. We replicate the scalar values.
    
    % Columns for diagonals: [Lower, Main, Upper]
    B = [val_sub * ones(K, 1), val_main * ones(K, 1), val_sup * ones(K, 1)];
    
    % Positions: -1 (sub), 0 (main), 1 (super)
    d_pos = [-1, 0, 1];
    
    % Create the sparse matrix
    P = spdiags(B, d_pos, K, K);
    
    % Note: If strictly cyclic preconditioner is needed (including corners),
    % one would add diagonals at positions -(K-1) and (K-1).
    % However, for efficient LU decomposition as requested, 
    % we strictly keep the central band.
end
