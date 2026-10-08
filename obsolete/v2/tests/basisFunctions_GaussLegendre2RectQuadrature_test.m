clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 3;
countY = 2;
sizeX = 2;
sizeY = 1;

rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
rect = rooftopMoM.mesh.UniformRect(rectPrimitive);
mesh = rect.generateMesh();

order = 3;

quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);
quad.generate()

quad.memory().bytes()
