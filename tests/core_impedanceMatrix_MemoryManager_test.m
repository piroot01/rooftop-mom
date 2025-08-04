clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 1;
countY = 100;
sizeX = 1;
sizeY = 1;

rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
rect = rooftopMoM.mesh.UniformRect(rectPrimitive);
mesh = rect.generateMesh();

order = 5;
k = 1;

quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);
rooftops = rooftopMoM.basisFunctions.UniformRooftops(mesh, quad);

rooftops.generate();

context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);

maxMemoryGB = 10;
maxMemory = rooftopMoM.Memory(maxMemoryGB, 'GB');

manager = rooftopMoM.core.impedanceMatrix.MemoryManager(maxMemory);
manager.registerContext(context);

elementCount = manager.elementCountMagneticPart(context);

planner = rooftopMoM.core.impedanceMatrix.AssemblerChunkPlanner(elementCount);

context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);

fullAssembler = rooftopMoM.core.impedanceMatrix.MagneticPartAssembler(planner);

Z = fullAssembler.generateFull(context);
