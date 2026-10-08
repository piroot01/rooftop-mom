
clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 3;
countY = 3;
sizeX = 1;
sizeY = 1;

rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);

order = 5;

quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);
rooftops = rooftopMoM.basisFunctions.FullRectUniformRooftops(rectPrimitive, quad);

rooftops.generate();

elementCount = 10000;

planner = rooftopMoM.core.greensMatrix.AssemblerChunkPlanner(elementCount);

k = 1;

context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);

chunksX = planner.generateDir(context, 'x');

celldisp(chunksX)
