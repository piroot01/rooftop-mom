clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 3;
countY = 3;
sizeX = 1;
sizeY = 1;

rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);

order = 2;

quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);

quad.generate();

points = rooftopMoM.basisFunctions.utils.generateQuadraturePointsUniformRectPrimitive(quad, rectPrimitive)
