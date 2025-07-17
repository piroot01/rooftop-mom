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

quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);
rooftops = rooftopMoM.basisFunctions.UniformRooftops(mesh, quad);

rooftops.generate();

elementCount = 1e7;

planner = rooftopMoM.core.impedanceMatrix.AssemblerChunkPlanner(elementCount);

lambda = 2 * countY * sizeY;
k = 2 * pi / lambda;

context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);

fullAssembler = rooftopMoM.core.impedanceMatrix.FullAssembler(planner);

Z = fullAssembler.generateFull(context);

mode = rooftopMoM.operators.utils.generateMiddlePortMode(rooftops);

V = 1;

excitation = rooftopMoM.operators.DeltaGapFeederExcitation.generate(rooftops, mode, V);

currents = rooftopMoM.core.solvers.Inversion.solve(Z, excitation);

Zin = rooftopMoM.operators.InputImpedance.generate(Z, currents, mode)
