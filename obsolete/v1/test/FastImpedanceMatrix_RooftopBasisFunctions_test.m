clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 2;
countY = 2;
gridSize = [countX, countY];

order = 2;

quadrature = mom.Quadrature(order);

rooftopBasisFunctions = mom.FastImpedanceMatrix.RooftopBasisFunctions(quadrature, gridSize);
rooftopBasisFunctions.generatePoints();

rooftopBasisFunctions.getRooftopPoints()

rooftopBasisFunctions.getQuadraturePoints()
