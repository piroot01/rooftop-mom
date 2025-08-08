clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 4;
countY = 3;
sizeX = 1.2;
sizeY = 0.8;

rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);

order = 2;

quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);

rooftops = rooftopMoM.basisFunctions.FullRectUniformRooftops(rectPrimitive, quad);

rooftops.generate();

mask = rooftopMoM.core.greensMatrix.utils.generateSingularityMask(rooftops.indicesDirection('x'))
