clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 3;
countY = 3;
sizeX = 2;
sizeY = 1;

rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY, genome);

rect = rooftopMoM.mesh.UniformRect(rectPrimitive);

mesh = rect.generateMesh();

mesh.faces()

rooftopMoM.mesh.utils.generateInteriorEdges(mesh.vertices(), mesh.faces())
