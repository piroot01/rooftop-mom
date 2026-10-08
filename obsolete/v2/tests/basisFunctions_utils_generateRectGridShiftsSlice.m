clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = uint32(3);
countY = uint32(3);
sizeX = 1.5;
sizeY = 0.5;

% full shifts
shifts = rooftopMoM.basisFunctions.utils.generateRectGridShifts(countX - 1, countY, sizeX, sizeY)

slice = 3;
shifts = rooftopMoM.basisFunctions.utils.generateRectGridShiftsSlice(countX - 1, countY, sizeX, sizeY, slice)
