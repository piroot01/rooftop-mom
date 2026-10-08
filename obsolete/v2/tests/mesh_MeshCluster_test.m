clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 3;
countY = 3;
sizeX = 2;
sizeY = 1;

rectPrimitive1 = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
rect1 = rooftopMoM.mesh.UniformRect(rectPrimitive1);
mesh1 = rect1.generateMesh();

origin = [10; 20; 30];

rectPrimitive2 = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
rect2 = rooftopMoM.mesh.UniformRect(rectPrimitive2, origin, 0, 0, 0);
mesh2 = rect2.generateMesh();

cluster = rooftopMoM.mesh.MeshCluster();
cluster.insert(mesh1);
cluster.insert(mesh2);

celldisp(cluster.vertices())
celldisp(cluster.faces())

disp(cluster.memory().bytes())
