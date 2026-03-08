clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

loggerDir = fullfile(currentDir, '..', '3rd', 'matlab-logger', 'src');
addpath(loggerDir);

% Configure the logger.
logger.Logger.configure( ...
    'timeFormat', 'ss.SSS', ...
    'minLevel', logger.LogLevel.TRACE, ...
    'sinks', { ...
        logger.sink.ConsoleSink(), ...
    } ...
);

countX = 3;
countY = 3;
sizeX = 1;
sizeY = 1;

order = 3;
elementCount = 1e6;
k = pi / (sizeY * countY);
V = 1;

quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);
rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);

rooftops = rooftopMoM.basisFunctions.UniformRectPrimitiveUniformRooftops(rectPrimitive, quad);
rooftops.generate();

planner = rooftopMoM.core.greensMatrix.AssemblerChunkPlanner(elementCount);
context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);
assembler = rooftopMoM.core.greensMatrix.FullAssembler(planner);

greensMatrix = assembler.generate(context);

spectralGreensMatrix = greensMatrix.generateSpectralGreensMatrix();

greensMatrixApplicator = rooftopMoM.core.greensMatrix.FFTAcceleratedApplicator( ...
    context, spectralGreensMatrix, greensMatrix.rooftopContent());

mode = rooftopMoM.operator.utils.generateMiddlePortMode(rooftops.count());
excitation = rooftopMoM.operator.DeltaGapFeederExcitation.generate(rooftops.count(), mode, V);

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
fromRooftops = toRooftops.';

% Now we tst if the matrix toRooftops is an orthogonal matrix.
fromRooftops * toRooftops

%{
% test of the hypothesis that H_{Lambda Sigma} from the paper is an transformation from Loop/Stars to Rooftops
e = zeros(rooftops.count(), 1);
e(1) = 1;

Lambda * e;
%}

tic;
LambdaInv = inv(Lambda);
toc;

% Tweak the number of start to match the count of the rooftops
rooftopCount = rooftops.count();
loopCount = size(loops, 2);
starCount = size(stars, 2);
removeCount = rooftopCount - loopCount - starCount;

S = starsTransformation(1:(end + removeCount), :).';

% now we construct the co-tree
coTreeBasisFunctions = setdiff(1:rooftopCount, abs(treeConverted));
coTreeSize = numel(coTreeBasisFunctions);

cols = 1:coTreeSize;
rows = coTreeBasisFunctions;
vals = ones(1, coTreeSize);

CT = sparse(rows, cols, vals, rooftopCount, loopCount);

% must be both zero
%CT.' * T
%S.' * L

tic;
LambdaInvFast = [(CT.' * L) \ CT.'; (S.' * T) \ S.'];
toc;

% Testing vector for the top line
topLine = (CT.' * L);

u = randi([-1, 1], size(topLine, 1), 1);

tic;
y1 = topLine \ u;
toc;

function y = fastTopLineMult(v, m, n)
    % fastBlockMultiply: Multiplies the structured block matrix by vector v.
    %
    % Inputs:
    %   v: Column vector of length m * n
    %   m: Size of each subblock (e.g., 4)
    %   n: Number of subblocks on the diagonal (e.g., 3)
    %
    % Output:
    %   y: Resulting vector
    
    y = zeros(size(v)); % Preallocate output
    
    for b = 1:n
        % Calculate global start and end indices for the current block
        startIdx = (b - 1) * m + 1;
        endIdx = b * m;
        
        if mod(b, 2) == 1
            % ODD BLOCKS: Lower triangular of -1s
            % Iterate forwards
            accumulator = 0;
            for i = startIdx:endIdx
                accumulator = accumulator - v(i);
                y(i) = accumulator;
            end
            
        else
            % EVEN BLOCKS: Upper triangular of 1s
            % Iterate backwards
            accumulator = 0;
            for i = endIdx:-1:startIdx
                accumulator = accumulator + v(i);
                y(i) = accumulator;
            end
        end
    end
end

tic;
y2 = fastTopLineMult(u, countX - 1, countY - 1);
toc;

all(y1 == y2)

% Testing vector for the bottom line
bottomLine = (S.' * T);

v = randi([-1, 1], size(bottomLine, 1), 1);

full(bottomLine)

tic;
x1 = bottomLine \ v;
toc;

function list = pivotList(tree)
    list = tree(2:(end - 1));
end

function x = applyBottomLineInv(list, b)
    N = numel(b);
    x = zeros(N, 1);

    x(1) = b(1);

    for c = 2:N
        r = list(c - 1);
        if r > N
            x(N - 1) = (-1) * (b(N) + b(N - 1));
            x(N) = (-1) * b(N - 1);
            break;
        else
            x(c) = x(c - 1) + b(r);
        end
    end
end

pivotList(tree)

tic;
x2 = applyBottomLineInv(pivotList(tree), v);
toc;

all(x1 == x2)

%{
T = full(treeTransformation).';

% Try to fix the orthogonality issue.
treeConvertedAbsSort = sort(abs(treeConverted));

rooftopCount = rooftops.count();

Talt = eye(rooftopCount, rooftopCount);
Talt = Talt(:, treeConvertedAbsSort);

coTree = setdiff(1:rooftopCount, abs(treeConverted));
Tco = eye(rooftopCount, rooftopCount);
Tco = Tco(:, coTree);

L = full(loopTransformation).';

% check if Tco^H * Talt == 0
fstCond = Tco' * T;
all(fstCond == 0, 'all')

S = full(starsTransformation).';

% to perfectly complement the loop basis functions we must drop some stars.
S = S(:, 1:(rooftopCount - size(L, 2)));

% check if S' * L == 0
sndCond = S' * L;
all(sndCond == 0, 'all')

% we see that both conditions are met, we can find the inverse of Lambda
Lambda = [L, T];

topLine = Tco' * L;
bottomLine = S' * T;

bottomLineSparse = sparse(bottomLine);

v = ones(size(bottomLine, 1), 1);

x1 = bottomLine \ v;

function list = pivotList(tree)
    list = tree(2:(end - 1));
end

function x = applyBottomLineInv(list, b)
    N = numel(b);
    x = zeros(N, 1);

    x(1) = b(1);

    for c = 2:N
        r = list(c - 1);
        if r > N
            x(N - 1) = (-1) * (b(N) + b(N - 1));
            x(N) = (-1) * b(N);
            break;
        else
            x(c) = x(c - 1) + b(r);
        end
    end
end

tic;
x2 = applyBottomLineInv(pivotList(tree), v);
toc;

tic;
x3 = bottomLineSparse \ v;
toc;

x1Andx2Same = all(x1 == x2)
x2Andx3Same = all(x1 == x3)
% now we can perform the inversion in O(N)

LambdaInv = [inv(topLine) * Tco'; inv(bottomLine) * S'];

inv(bottomLine);

all(inv(Lambda) == LambdaInv, 'all');
%}

%excitationLoopTree = Lambda.' * excitation;
excitationLoopTree = excitation;

function y = apply(x, greensMatrixApplicator, T, Tinv)
    y = Tinv * greensMatrixApplicator.apply(T * x);
    %y = greensMatrixApplicator.apply(x);
end

applyHandle = @(x) apply(x, greensMatrixApplicator, Lambda, Lambda.');

solver = rooftopMoM.core.solver.ConjugateGradient(applyHandle);

tolerance = 1e-6;
maxIterCount = 10000;

outputConvSolver = solver.solve( ...
    'excitation', excitationLoopTree, ...
    'tolerance', tolerance, ...
    'maxIterCount', maxIterCount ...
);

fprintf('Number of the iterations: %d\n', outputConvSolver.count());
