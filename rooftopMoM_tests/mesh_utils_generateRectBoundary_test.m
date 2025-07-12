clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 3;
countY = 3;
sizeX = 2;
sizeY = 1;

count = countX * countY;

missing = randperm(count, 4);

genome = true(count, 1);
genome(missing) = false;

rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);

rect = rooftopMoM.mesh.UniformRect(rectPrimitive);

mesh = rect.generateMesh();

rooftopMoM.mesh.utils.generateRectBoundary(mesh.vertices(), mesh.faces())
