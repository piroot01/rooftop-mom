clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 3;
countY = 3;
sizeX = 2;
sizeY = 1;

rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);

origin = [1; 2; 3];

rect = rooftopMoM.mesh.UniformRect(rectPrimitive, origin, pi/3, 0, 0);

mesh = rect.generateMesh();

[s1, s2] = rooftopMoM.mesh.utils.generateRectSize(mesh.vertices{1}, mesh.faces{1})
