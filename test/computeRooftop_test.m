clc; clearvars;

currentFolder = fileparts(mfilename('fullpath'));
addpath(currentFolder, '..', 'src/');

countX = 3;
countY = 3;
sizeX = 0.5;
sizeY = 1.5;
origin = [0, 0, 2];

rectPrim = mom.MeshRectPrimitive(countX, countY, sizeX, sizeY, []);
rectMesh = mom.MeshRect(rectPrim, origin, 0, 0, 0);

mesh = rectMesh.getAsMesh();

basisFunctions = mom.BasisFunctions(mesh);
basisFunc = basisFunctions.get();

order = 2;
order2 = order^2;

basisFuncCount = size(basisFunc, 2);
psiP = zeros(3, order^2 * basisFuncCount);
psiN = zeros(3, order^2 * basisFuncCount);

basisFuncIdx = 1:basisFuncCount;

quadPoints = mom.QuadPoints(order);

for n = basisFuncIdx
    chosenOne = basisFunc(:, n);
    vertexIndices = mesh.faces(:, chosenOne);
    face = mesh.vertices(:, vertexIndices(:));
    qP = quadPoints.getTransformed(face);
    [rP, rN] = mom.computeRooftop(face, qP);
    psiP(:, (order2 * (n - 1) + 1):(order2 * n)) = rP;
    psiN(:, (order2 * (n - 1) + 1):(order2 * n)) = rN;
end

disp(psiP);
disp(psiN);