clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 3;
countY = 3;
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

accelerator = rooftopMoM.core.impedanceMatrix.InteractionSymmetriesAccelerator(fullAssembler);

Z = accelerator.generateFull(context);
