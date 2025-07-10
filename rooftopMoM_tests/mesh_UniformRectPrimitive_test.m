clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 3;
countY = 3;
sizeX = 2;
sizeY = 1;

missing = randi([1, countX * countY]);

genome = true(countX * countY, 1);
genome(missing) = false;

rect = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY, genome);

rect.indices()
rect.indicesGenome()
