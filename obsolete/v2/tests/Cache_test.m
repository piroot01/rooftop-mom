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

% Enable caching
rooftopMoM.Cache.setStatus(true);
rooftopMoM.Cache.setOption(rooftopMoM.CacheKeys.CORE_IMPEDANCEMATRIX_FULLASSEMBLER_ROOFTOPPOINTS, true);

countX = 3;
countY = 3;
sizeX = 1;
sizeY = 1;

rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
rect = rooftopMoM.mesh.UniformRect(rectPrimitive);
mesh = rect.generateMesh();

order = 2;

quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);
rooftops = rooftopMoM.basisFunctions.MeshUniformRooftops(mesh, quad);

rooftops.generate();

elementCount = 100000;

planner = rooftopMoM.core.impedanceMatrix.AssemblerChunkPlanner(elementCount);

k = 1;

context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);

fullAssembler = rooftopMoM.core.impedanceMatrix.FullAssembler(planner);

Z = fullAssembler.generateFull(context);

rooftopMoM.Cache.getOptionStatus(rooftopMoM.CacheKeys.CORE_IMPEDANCEMATRIX_FULLASSEMBLER_ROOFTOPPOINTS)

% Obtain the rooftopPoints
rooftopPoints = rooftopMoM.Cache.get(rooftopMoM.CacheKeys.CORE_IMPEDANCEMATRIX_FULLASSEMBLER_ROOFTOPPOINTS);

pre = rooftopMoM.Cache.getOptionDataStatus(rooftopMoM.CacheKeys.CORE_IMPEDANCEMATRIX_FULLASSEMBLER_ROOFTOPPOINTS)

rooftopMoM.Cache.clear(rooftopMoM.CacheKeys.CORE_IMPEDANCEMATRIX_FULLASSEMBLER_ROOFTOPPOINTS);

rooftopMoM.Cache.getOptionStatus(rooftopMoM.CacheKeys.CORE_IMPEDANCEMATRIX_FULLASSEMBLER_ROOFTOPPOINTS)

post = rooftopMoM.Cache.getOptionDataStatus(rooftopMoM.CacheKeys.CORE_IMPEDANCEMATRIX_FULLASSEMBLER_ROOFTOPPOINTS)

% Destory the cache
rooftopMoM.Cache.destroy();

rooftopMoM.Cache.getOptionStatus(rooftopMoM.CacheKeys.CORE_IMPEDANCEMATRIX_FULLASSEMBLER_ROOFTOPPOINTS)
