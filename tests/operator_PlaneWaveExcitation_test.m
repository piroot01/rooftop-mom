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

countX = 5000;
countY = 5000;
sizeX = 1;
sizeY = 1;

k = 1;
order = 2;
elementCount = 1e7;

direction = [0; 0; -1];
polarization = [0; 1; 0];

quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);

rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
rooftops = rooftopMoM.basisFunctions.UniformRectPrimitiveUniformRooftops(rectPrimitive, quad);
rooftops.generate();
context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);

tic;
V2 = rooftopMoM.operator.PlaneWaveExcitation.generateWithAssemblerContext(context, direction, polarization, 1, elementCount);
toc;
