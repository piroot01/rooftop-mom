clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

loggerDir = fullfile(currentDir, '..', '3rd', 'matlab-logger', 'src');
addpath(loggerDir);

% Configure the logger.
logger.Logger.configure( ...
    'timeFormat', 'ss.SSS', ...
    'minLevel', logger.LogLevel.DEBUG, ...
    'sinks', { ...
        logger.sink.ConsoleSink(), ...
    } ...
);

% Setup caching
rooftopMoM.Cache.setStatus(true);
rooftopMoM.Cache.setCacheAll(true);

countX = 100;
countY = 100;
sizeX = 1.5;
sizeY = 0.5;

k = 1;
order = 5;
elementCount = 1e7;

quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);
%{
%% Impedance matrix
rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
rect = rooftopMoM.mesh.UniformRect(rectPrimitive);
mesh = rect.generateMesh();

rooftops = rooftopMoM.basisFunctions.MeshUniformRooftops(mesh, quad);

rooftops.generate();

planner = rooftopMoM.core.impedanceMatrix.AssemblerChunkPlanner(elementCount);

context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);

fullAssembler = rooftopMoM.core.impedanceMatrix.FullAssembler(planner);

Z = fullAssembler.generateFull(context);

tic;
[qP, rP, rV] = rooftopMoM.operator.utils.generateRooftopData(rooftops, true);
toc;

rooftopMoM.Cache.clear();

tic;
[qP, rP, rV] = rooftopMoM.operator.utils.generateRooftopData(rooftops, false);
toc;
%}

%% Greens matrix
rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
rooftops = rooftopMoM.basisFunctions.UniformRectPrimitiveUniformRooftops(rectPrimitive, quad);

rooftops.generate();

planner = rooftopMoM.core.greensMatrix.AssemblerChunkPlanner(elementCount);

context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);

assembler = rooftopMoM.core.greensMatrix.FullAssembler(planner);

output = assembler.generateFull(context);

tic;
[qP, rP, rV] = rooftopMoM.operator.utils.generateRooftopData(rooftops, true);
toc;

rooftopMoM.Cache.clear();

tic;
[qP, rP, rV] = rooftopMoM.operator.utils.generateRooftopData(rooftops, false);
toc;
