clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

function result = runTest(countX, countY, sizeX, sizeY)
    rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
    rect = rooftopMoM.mesh.UniformRect(rectPrimitive);
    mesh = rect.generateMesh();

    order = 5;

    quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);
    rooftops = rooftopMoM.basisFunctions.UniformRooftops(mesh, quad);

    rooftops.generate();

    elementCount = 1e6;

    planner = rooftopMoM.core.impedanceMatrix.AssemblerChunkPlanner(elementCount);

    k = 1;

    context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);

    fullAssembler = rooftopMoM.core.impedanceMatrix.FullAssembler(planner);

    Z_non = fullAssembler.generateFull(context);

    % Compute using the accelerator.
    accelerator = rooftopMoM.core.impedanceMatrix.InteractionSymmetriesAccelerator(fullAssembler);

    Z_acc = accelerator.generateFull(context);

    result = max(Z_non - Z_acc, [], 'all');
end

% 1. test: Square radiator square grid.
test1Verdict = runTest(3, 3, 1, 1)

% 2. test: Rect Y longer square grid.
test2Verdict = runTest(3, 4, 1, 1)

% 3. test: Rect X longer square grid.
test3Verdict = runTest(4, 3, 1, 1)

% 4. test: Square radiator rect grid.
test4Verdict = runTest(3, 3, 1.2, 1)

% 5. test: Rect Y longer rect grid.
test5Verdict = runTest(3, 4, 1.2, 1)

% 6. test: Rect X longer rect grid.
test6Verdict = runTest(4, 3, 1.2, 1)

% 7. test: Dipole X
test7Verdict = runTest(10, 1, 1, 1)

% 8. test: Dipole Y
test8Verdict = runTest(1, 10, 1, 1)
