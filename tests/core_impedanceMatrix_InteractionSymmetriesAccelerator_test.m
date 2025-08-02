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

% 1. test: Square qrid square radiator.
test1Verdict = runTest(2, 3, 1, 1)
