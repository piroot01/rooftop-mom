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

% Now we add the basis function back.
reorderIndices = shiftAndInsert(allIndices_, selectedBasisFunction);

Yuu = rank1Add(Yu_cut, Z, allIndices_, reorderIndices, selectedBasisFunction);

if all(Yuu - Y < 100 * eps, 'all')
    fprintf("The matrices Yuu and Y are identical.\n");
end

% remove another basis function.
another = 6;
Yu = rank1Removal(Yu, another);

% and now we will bring back both.
Yuu_zero = rank1Add_zero(Yu, Z, [selectedBasisFunction, another], selectedBasisFunction);
Yuu_zero = rank1Add_zero(Yuu_zero, Z, another, another);

if all(Yuu_zero - Yuu < 100 * eps, 'all')
    fprintf("The matrices Yuu_zero and Yuu are identical.\n");
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

function Y = rank1Add_zero(Y, Z, zeroIdx, n)
% Y         (N×N) current inverse with zeros on rows/cols in zeroIdx
% Z         (N×N) original matrix
% zeroIdx   list of indices where Y is currently empty
% n         the specific index in zeroIdx to “fill in” now
%
% After this call, row n and col n of Y will have been updated
% by the Woodbury rank‑1 formula, and n is removed from zeroIdx.

    % build the list of the “old” active indices
    activeIdx = setdiff(1:size(Y,1), zeroIdx);

    %--- 1) compute scalar and cross‐vector from Z
    Znn = Z(n,n);
    zn  = Z(activeIdx, n);               % (m×1) cross‐columns

    %--- 2) x_n = Y(active,active) * z_n
    xn  = Y(activeIdx, activeIdx) * zn;  % (m×1)

    %--- 3) compute Woodbury factor
    cn  = 1 / (Znn - zn.' * xn);         % scalar

    %--- 4) update the active‐active block in place
    Y(activeIdx, activeIdx) = ...
        Y(activeIdx, activeIdx) + cn * (xn * xn.');

    %--- 5) fill in the new off‐diagonals
    Y(activeIdx, n) = -cn * xn;
    Y(n, activeIdx) =  Y(activeIdx, n).';

    %--- 6) set the new diagonal
    Y(n, n) =  cn;

    % (Optionally) if you need to keep track of which indices remain empty,
    % you’d remove n from zeroIdx here, e.g.
    % zeroIdx(zeroIdx==n) = [];
end
