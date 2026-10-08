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

countX = 5000;
countY = 2500;
sizeX = 1;
sizeY = 1;

k = 1;
order = 3;
elementCount = 3e7;

direction = [0; 0; -1];
polarization = [0; 1.5; 0];

quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);

rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
rooftops = rooftopMoM.basisFunctions.UniformRectPrimitiveUniformRooftops(rectPrimitive, quad);
rooftops.generate();
context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);

tic;
V1 = rooftopMoM.operator.PlaneWaveExcitation.generateWithAssemblerContext(context, direction, polarization, 1, elementCount);
toc;

tic;
V2 = rooftopMoM.operator.PlaneWaveExcitation.generatePerpendicular(rooftops, polarization, 1);
toc;

all(V1 - V2 < 1e-10)