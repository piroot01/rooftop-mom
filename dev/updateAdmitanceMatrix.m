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

order = 5;

quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);
rooftops = rooftopMoM.basisFunctions.UniformRooftops(mesh, quad);

rooftops.generate();

elementCount = 100000;

planner = rooftopMoM.core.impedanceMatrix.AssemblerChunkPlanner(elementCount);

k = 1;

context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);

fullAssembler = rooftopMoM.core.impedanceMatrix.FullAssembler(planner);

Z = fullAssembler.generateFull(context);

Y = inv(Z);

pertAdmittanceMatrix = PerturbedMatrix(Y)

N     = 100000000;
a     = sort(randi(1e6, 1, N));
k     = randi([0, N], 1, 1);

tic;

[~, retrace] = sort([a, k]);

toc;

tic;

fastestRetrace = shiftAndInsert(a, k);

toc;

all(retrace == fastestRetrace)

function perturbedMatrix = rank1Removal(matrix, index)
    matrixColumn = matrix(:, index);
    matrixDiagElement = matrix(index, index);

    perturbedMatrix = matrix - (matrixColumn / matrixDiagElement) * matrixColumn.';
end

function Y_updated = rank1Add(Y, Z_original, indices, n)

    Znn = Z_original(n, n);
    zn = Z_original(indices, n);
    xn = Y * zn;
    cn = 1 / (Znn - zn.' * xn);
    cn_xn = cn .* xn;
    Y_updated = [ Y + cn_xn * xn.', (-1) * cn_xn; ...
                  (-1) * xn_xn.', cn ];
end
