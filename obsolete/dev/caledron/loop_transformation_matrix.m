clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', '..', 'src');
addpath(srcDir);

loggerDir = fullfile(currentDir, '..', '..', '3rd', 'matlab-logger', 'src');
addpath(loggerDir);

% Configure the logger.
logger.Logger.configure( ...
    'timeFormat', 'ss.SSS', ...
    'minLevel', logger.LogLevel.TRACE, ...
    'sinks', { ...
        logger.sink.ConsoleSink(), ...
    } ...
);

countX = 3;
countY = 3;
sizeX = 1;
sizeY = 1;

order = 3;
elementCount = 1e6;
k = pi / (sizeY * countY);
V = 1;

quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);
rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);

basisFunctionsIndices = rooftopMoM.basisFunctions.UniformRectPrimitiveDivergenceConforming(rectPrimitive);
basisFunctionsIndices.generate();
indices = basisFunctionsIndices.indices()

% total basis function count
N = basisFunctionsIndices.count();
Nx = basisFunctionsIndices.countDir('x');
Ny = basisFunctionsIndices.countDir('y');

% total count of the loops
% TODO: get from the N
Nl = (countX - 1) * (countY - 1);

% first we define an empty L matrix of size Nl x N
L = zeros(Nl, N);

% fill the matrix
for cornerIndex = 1:Nl
    loopIndices = [cornerIndex, cornerIndex + Nx + 1, cornerIndex + countX - 1, cornerIndex + Nx]

    L(cornerIndex, loopIndices(1:2)) = 1;
    L(cornerIndex, loopIndices(3:4)) = -1;
end

L