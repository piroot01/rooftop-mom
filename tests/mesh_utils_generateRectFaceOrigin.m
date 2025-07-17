clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 3;
countY = 3;
sizeX = 2;
sizeY = 1;

rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
rect = rooftopMoM.mesh.UniformRect(rectPrimitive);
mesh = rect.generateMesh();

selectedFaces = [1; 5; 9];

rooftopMoM.mesh.utils.generateRectFaceOrigin(mesh.vertices(), mesh.faces(), 4)

rooftopMoM.mesh.utils.generateRectFaceOrigin(mesh.vertices(), mesh.faces(), selectedFaces)
