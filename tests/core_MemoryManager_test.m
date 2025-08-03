clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

memoryMB = 10;

manager = rooftopMoM.core.MemoryManager(memoryMB, 'MB');
