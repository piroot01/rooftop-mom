clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 3;
countY = 2;
sizeX = 0.5;
sizeY = 1;

rectPrimitive1 = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
rect1 = rooftopMoM.mesh.UniformRect(rectPrimitive1);
mesh1 = rect1.generateMesh();

origin = [10; 20; 30];

rectPrimitive2 = rooftopMoM.mesh.UniformRectPrimitive(countX * 2, countY, sizeX, sizeY);
rect2 = rooftopMoM.mesh.UniformRect(rectPrimitive2, origin, 0, 0, 0);
mesh2 = rect2.generateMesh();

cluster = rooftopMoM.mesh.MeshCluster();
cluster.insert(mesh1);
cluster.insert(mesh2);

order = 5;

quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);
rooftops = rooftopMoM.basisFunctions.UniformRooftops(cluster, quad);

rooftops.generate();
[points, values] = rooftops.generateRooftops()

points{1}

rooftops.memory().bytes()

mem = rooftopMoM.basisFunctions.utils.memoryForRooftopPoints(rooftops);

mem.bytes()

% genome - sloppy but should work
%{
genome = true(1, rectPrimitive2.count());
genome(3) = false;

mesh2Genome = rooftopMoM.mesh.utils.filterMeshByGenome(mesh2, genome);

cluster.remove(2);
cluster.insert(mesh2Genome)

rooftopsGenome = rooftopMoM.basisFunctions.UniformRooftops(cluster, quad);

rooftops.generate();
rooftops.generateQuadraturePoints()
%}
