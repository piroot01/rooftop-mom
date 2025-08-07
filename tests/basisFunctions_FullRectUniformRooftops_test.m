clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 1;
countY = 2;
sizeX = 1.2;
sizeY = 0.8;

rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
rect = rooftopMoM.mesh.UniformRect(rectPrimitive, [0; 0; 0], 0, 0, 0);
mesh = rect.generateMesh();

order = 2;

quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);

rooftops = rooftopMoM.basisFunctions.UniformRooftops(mesh, quad);
rooftops.generate();
rooftops.generateRooftops();

rooftopsFull = rooftopMoM.basisFunctions.FullRectUniformRooftops(rectPrimitive, quad);

rooftopsFull.generate();

rooftopsFull.generateRooftops();
