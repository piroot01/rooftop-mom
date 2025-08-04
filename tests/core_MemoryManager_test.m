clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

memoryMB = 10;

maxMem = rooftopMoM.Memory(memoryMB, 'MB');

order = 4;

quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);

manager = rooftopMoM.core.MemoryManager(maxMem);

manager.register(quad);
