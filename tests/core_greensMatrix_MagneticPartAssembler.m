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
rooftops = rooftopMoM.basisFunctions.UniformRectPrimitiveUniformRooftops(rectPrimitive, quad);

rooftops.generate();

elementCount = 100 * order^4;

planner = rooftopMoM.core.greensMatrix.AssemblerChunkPlanner(elementCount);

k = 1;

context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);

magAssembler = rooftopMoM.core.greensMatrix.MagneticPartAssembler(planner);

%Gxx = magAssembler.generateFullDirection(context, 'x')
Gyy = magAssembler.generateFullDir(context, 'y')

% Test comparing to the impedance matrix
%{
rect = rooftopMoM.mesh.UniformRect(rectPrimitive);
mesh = rect.generateMesh();

rooftops = rooftopMoM.basisFunctions.UniformRooftops(mesh, quad);
rooftops.generate();

context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);

planner = rooftopMoM.core.impedanceMatrix.AssemblerChunkPlanner(elementCount);

magAssembler = rooftopMoM.core.impedanceMatrix.MagneticPartAssembler(planner);

Z = magAssembler.generateFull(context);
%}
