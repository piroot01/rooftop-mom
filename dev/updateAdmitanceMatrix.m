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

% test setup: we remove and add selected basis function
selectedBasisFunction = 4;

% Determine the current active indices
allIndices = 1:rooftops.count();

% remove basis function
Yu = rank1Removal(Y, selectedBasisFunction);
allIndices(selectedBasisFunction) = [];

Yu = Yu(allIndices, allIndices);

% Check if the Yu is correct
genome{1} = true(size(Z, 2), 1);
genome{1}(selectedBasisFunction) = false;

rooftops.generate(genome);

Zu = fullAssembler.generateFull(context);

Yu_ = inv(Zu);

if all(Yu - Yu_ < 10 * eps, 'all')
    fprintf("The matrices Yu and Yu_ are identical.\n");
end

% Now we add the basis function back.
reorderIndices = shiftAndInsert(allIndices, selectedBasisFunction);

Yuu = rank1Add(Yu, Z, allIndices, reorderIndices, selectedBasisFunction)

if all(Yuu - Y < 10 * eps, 'all')
    fprintf("The matrices Yuu and Y are identical.\n");
end

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
