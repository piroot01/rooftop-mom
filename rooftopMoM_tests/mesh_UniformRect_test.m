clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 3;
countY = 3;
sizeX = 2;
sizeY = 1;

missing = randperm(countX * countY, 6)

genome = true(countX * countY, 1);
genome(missing) = false;

rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY, genome);

rect = rooftopMoM.mesh.UniformRect(rectPrimitive);

rect.generateMesh();
