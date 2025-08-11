clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 4;
countY = 4;
sizeX = 1;
sizeY = 1;

rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);

order = 5;

quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);
rooftops = rooftopMoM.basisFunctions.UniformRectPrimitiveUniformRooftops(rectPrimitive, quad);

rooftops.generate();

elementCount = 100 * order^4;

planner = rooftopMoM.core.greensMatrix.AssemblerChunkPlanner(elementCount);

k = 1;

context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);

elAssembler = rooftopMoM.core.greensMatrix.ElectricPartAssembler(planner);

Ge = elAssembler.generateFull(context)
