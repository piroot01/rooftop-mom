clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 1;
countY = 10;
sizeX = 1;
sizeY = 1;

rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
rect = rooftopMoM.mesh.UniformRect(rectPrimitive);
mesh = rect.generateMesh();

order = 10;

quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);
rooftops = rooftopMoM.basisFunctions.UniformRooftops(mesh, quad);

rooftops.generate();

elementCount = 1e7;

planner = rooftopMoM.core.impedanceMatrix.AssemblerChunkPlanner(elementCount);

k = 1;

context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);

fullAssembler = rooftopMoM.core.impedanceMatrix.FullAssembler(planner);

Z = fullAssembler.generateFull(context);

Y = inv(Z);

% test setup: we remove and add selected basis function
selectedBasisFunction = 4;

% Determine the current active indices
allIndices = 1:rooftops.count();

% remove basis function
Yu = rank1Removal(Y, selectedBasisFunction);
allIndices_ = allIndices;
allIndices_(selectedBasisFunction) = [];

Yu_cut = Yu(allIndices_, allIndices_);

% Check if the Yu is correct
genome{1} = true(size(Z, 2), 1);
genome{1}(selectedBasisFunction) = false;

rooftops.generate(genome);

Zu = fullAssembler.generateFull(context);

Yu_ = inv(Zu);

if all(Yu_cut - Yu_ < 10 * eps, 'all')
    fprintf("The matrices Yu_cut and Yu_ are identical.\n");
end

tic;
% Now we add the basis function back.
reorderIndices = shiftAndInsert(allIndices_, selectedBasisFunction);

Yuu = rank1Add(Yu_cut, Z, allIndices_, reorderIndices, selectedBasisFunction);
classical = toc

if all(Yuu - Y < 100 * eps, 'all')
    fprintf("The matrices Yuu and Y are identical.\n");
end

% remove another basis function.
another = 6;
Yu = rank1Removal(Yu, another);

activeIdx = 1:(countY - 1);
activeIdx([selectedBasisFunction, another]) = [];

% and now we will bring back both.
tic;
[Yuu_zero, activeIdx] = rank1AdditionInverseMatrix(Yu, Z, activeIdx, selectedBasisFunction);
newVersion = toc

[Yuu_zero, activeIdx] = rank1AdditionInverseMatrix(Yuu_zero, Z, activeIdx, another);

if all(Yuu_zero - Y < 100 * eps, 'all')
    fprintf("The matrices Yuu_zero and Y are identical.\n");
end

% now let us test the speed of each rank1Add
tries = 40;

count = 5000;
M = diag(rand(count, 1));
Minv = inv(M);
bf = 1;

allIdx = 1:count;

% the original variant
tic;

for n = 1:tries

MinvR = rank1Removal(Minv, bf);
allIdx(bf) = [];
MinvR = MinvR(allIdx, allIdx);

reorderIndices = shiftAndInsert(allIdx, bf);

MinvRA = rank1Add(MinvR, M, allIdx, reorderIndices, bf);

end

original = toc

allIdx = 1:count;

% new version

tic;

for n = 1:tries

[MinvR, allIdx] = rank1RemovalInverseMatrix(Minv, allIdx, bf);
[MinvRA, allIdx] = rank1AdditionInverseMatrix(MinvR, M, allIdx, bf);

end

newVersion = toc

function Y_updated = rank1Removal(Y, n)
    yn = Y(:, n);
    Ynn = Y(n, n);

    Y_updated = Y - (yn / Ynn) * yn.';
end

function Y_updated = rank1Add(Y, Z_original, indices, reorderIndices, n)
    Znn = Z_original(n, n);
    zn = Z_original(indices, n);
    xn = Y * zn;
    cn = 1 / (Znn - zn.' * xn);
    cn_xn = cn .* xn;
    Y_updated = [ Y + cn_xn * xn.', (-1) * cn_xn; ...
                  (-1) * cn_xn.', cn ];
    Y_updated = Y_updated(reorderIndices, reorderIndices);
end
