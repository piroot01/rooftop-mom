clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

loggerDir = fullfile(currentDir, '..', '3rd', 'matlab-logger', 'src');
addpath(loggerDir);

% Configure the logger.
logger.Logger.configure( ...
    'timeFormat', 'ss.SSS', ...
    'minLevel', logger.LogLevel.DEBUG, ...
    'sinks', { ...
        logger.sink.ConsoleSink(), ...
    } ...
);

countX = 100;
countY = 100;
sizeX = 1;
sizeY = 1;

order = 3;
elementCount = 1e6;
k = pi / (1000 * sizeY * countY);
%k = 1e-6;
V = 1;

quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);
rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);

rooftops = rooftopMoM.basisFunctions.UniformRectPrimitiveUniformRooftops(rectPrimitive, quad);
rooftops.generate();
%{
planner = rooftopMoM.core.greensMatrix.AssemblerChunkPlanner(elementCount);
context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);
assembler = rooftopMoM.core.greensMatrix.FullAssembler(planner);

greensMatrix = assembler.generate(context);

spectralGreensMatrix = greensMatrix.generateSpectralGreensMatrix();

greensMatrixApplicator = rooftopMoM.core.greensMatrix.FFTAcceleratedApplicator( ...
    context, spectralGreensMatrix, greensMatrix.rooftopContent());

mode = rooftopMoM.operator.utils.generateMiddlePortMode(rooftops.count());
excitation = rooftopMoM.operator.DeltaGapFeederExcitation.generate(rooftops.count(), mode, V);

% Reference solution using only the conjugate gradient method
solver = rooftopMoM.core.solver.ConjugateGradient(greensMatrixApplicator);

tolerance = 1e-12;
maxIterCount = 10000;

outputConvSolverRef = solver.solve( ...
    'excitation', excitation, ...
    'tolerance', tolerance, ...
    'maxIterCount', maxIterCount ...
);

currentRef = outputConvSolverRef.postRunInfo().current;
%}
% Prepare the loop-tree transformation matrix.
loops = rooftopMoM.core.solver.altBasis.generateLoops(countX, countY);
loopsConverted = rooftopMoM.core.solver.altBasis.convertLoopToRooftopIndices(loops, countX, countY);
loopTransformation = rooftopMoM.core.solver.altBasis.generateRooftopToLoopTransformation(loopsConverted, sizeX, sizeY);

tree = rooftopMoM.core.solver.altBasis.generateTree(countX, countY);
treeConverted = rooftopMoM.core.solver.altBasis.convertTreeToRooftopIndices(tree, countX, countY);
treeTransformation = rooftopMoM.core.solver.altBasis.generateRooftopToTreeTransformation(treeConverted, rooftops.count());

L = loopTransformation.';
T = treeTransformation.';

% Construct the coordinate transformation matrices.
toRooftops = cat(2, L, T);

% Frequency normalization
loopCount = size(L, 2);

function omega = toOmega(k)
    omega = k * physconst('LightSpeed');
end

omega0 = toOmega(k);

function y = applyPpr(x, omega, loopCount)
    y = x;
    y(1:loopCount) = y(1:loopCount) * (1 / omega);
end

function y = applyPpo(x, omega, loopCount)
    y = x;
    y((loopCount + 1):end) = y((loopCount + 1):end) * omega;
end

% basis rearragement
% 1. We assemble the K matrix K_{mn} = <P_m, div(R^{(T)}_N)>, where R^{(T)} are the tree basis functions.
treeSize = size(T, 2)

K = zeros(treeSize, treeSize);

for m = 1:treeSize
    for n = 1:treeSize
        % First we pick the two patches associated with the current tree rooftop.
        % The first patch will always be the positive one and the second one the negative one.
        patchIndices = tree(n:(n+1));

        % Since we are basically taking only one part of the rooftop based on the index of the
        % pulse function we can store just the polarity in the K matrix.
        % TODO: Add a proper rooftop divergence.
        if m == patchIndices(1)
            K(m, n) = 1;
        elseif m == patchIndices(2)
            K(m, n) = -1;
        end
    end
end

%disp(K);

% Now the idea is to be able to apply the inverse of such matrix in O(N) time and O(1) memory.
list = tree(2:(end - 1));
Q = randi([-1, 1], [treeSize, 1]);
tic;
I1 = K \ Q;
toc;

tic;
I2 = solveCurrentByList(list, Q, -1i);
toc;

all(I1 == I2)

%{
function I = solveCurrentByList(list, Q, omega)
    N = numel(Q);
    I = zeros(N, 1);
    I(1) = Q(1);

    lastRowHasMissing = false;
    lastRowHasMissingIdx = 0;

    for n = 2:N
        p = list(n - 1);
        if p > N
            lastRowHasMissing = true;
            lastRowHasMissingIdx = n;
            break;
        end
        I(n) = 1i * omega * Q(p) + I(n - 1);
    end

    if lastRowHasMissing
        I(N) = (-1i) * omega * Q(lastRowHasMissingIdx);
        k = 1;
        for n = (lastRowHasMissingIdx+1):N
            I(N - k) = I(N - k + 1) - 1i * omega * Q(n);
            k = k + 1;
        end
    end
end
%}

function I = solveCurrentByList(list, Q, omega)
    N = numel(Q);
    I = zeros(N, 1);
    
    if N == 0
        return;
    end
    
    I(1) = Q(1);
    
    if N == 1
        return; 
    end

    % Precalculate constant to avoid repeated complex multiplication
    i_omega = 1i * omega;

    % Find the first occurrence where a list value exceeds N.
    % We only search up to N-1, matching the original loop's scope.
    idx = find(list(1:N-1) > N, 1);

    if isempty(idx)
        % No missing row: purely a forward cumulative sum.
        list_slice = list(1:N-1);
        
        % Factoring i_omega out of the cumsum saves N-1 complex multiplications
        I(2:N) = Q(1) + i_omega * cumsum(Q(list_slice(:)));
    else
        break_n = idx + 1; % Matches original 'lastRowHasMissingIdx'

        % 1. Forward pass (computes up to break_n - 1)
        if break_n > 2
            list_slice = list(1:break_n-2);
            I(2:break_n-1) = Q(1) + i_omega * cumsum(Q(list_slice(:)));
        end

        % 2. Backward pass (starts from N down to break_n)
        I_N = -i_omega * Q(break_n);
        I(N) = I_N;

        if N > break_n
            % The original backward loop subtracts successive Q(n) values.
            % This is mathematically identical to subtracting a cumulative sum.
            Q_back = Q(break_n+1:N);
            I(N-1:-1:break_n) = I_N - i_omega * cumsum(Q_back(:));
        end
    end
end

function I = solveCurrent(Q, omega)
    N = numel(Q);
    I = zeros(N, 1);
    I(1) = Q(1);
    for n = 2:N
        I(n) = Q(n - 1) + 1i * omega * Q(n);
    end
end

function I = solveCurrentTranspose(Q)
    N = length(Q);
    I = zeros(N, 1);
    sum = 0;
    for n = N:(-1):1
        sum = sum + O(n);
        I(n) = sum;
    end
end

function y = apply(x, greensMatrixApplicator, toRooftops, omega, loopCount)
    % 1. Scale the normalized loop-tree trial vector (Apply P_PO)
    x_scaled = applyPpo(x, omega, loopCount);
    
    % 2. Transform from loop-tree to rooftop coefficients (Apply F)
    I_rt = toRooftops * x_scaled;
    
    % 3. Apply standard MoM Green's matrix (Apply Z_RT)
    V_rt = greensMatrixApplicator.apply(I_rt);
    
    % 4. Transform testing functions from rooftops to loop-trees (Apply F^T)
    V_lt = toRooftops.' * V_rt;
    
    % 5. Scale the resulting loop-tree testing vector (Apply P_PR)
    y = applyPpr(V_lt, omega, loopCount);
end

applyHandle = @(x) apply(x, greensMatrixApplicator, toRooftops, omega0, loopCount);

solver = rooftopMoM.core.solver.ConjugateGradient(applyHandle);
%{
excitationLoopTree = applyPpr(toRooftops.' * excitation, omega0, loopCount);

outputConvSolverLoopTree = solver.solve( ...
    'excitation', excitationLoopTree, ...
    'tolerance', tolerance, ...
    'maxIterCount', maxIterCount ...
);

currentLoopTree = outputConvSolverLoopTree.postRunInfo().current;

outputConvSolverRef.count()
outputConvSolverLoopTree.count()

max(currentRef - toRooftops * applyPpo(currentLoopTree, omega0, loopCount))
%}
