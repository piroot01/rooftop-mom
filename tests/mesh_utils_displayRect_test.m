clc; clearvars; close all;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 50;
countY = 50;
sizeX = 1;
sizeY = 1;

count = countX * countY;

missing = randperm(count, countX * countY / 2);

genome = true(count, 1);
genome(missing) = false;

rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);

rect = rooftopMoM.mesh.UniformRect(rectPrimitive);

mesh = rect.generateMesh();

rooftopMoM.mesh.utils.displayRect(mesh.vertices{1}, mesh.faces{1}, genome);
