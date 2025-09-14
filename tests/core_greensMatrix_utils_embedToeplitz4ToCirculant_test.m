clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

%
% test 1: cY == 1 and cX > 1, no padding
%

cX = randi(10, 1);
cY = 1;

Tgen1 = randi(100, [cX, cY])
H1 = rooftopMoM.core.greensMatrix.utils.embedToeplitz4ToCirculant(Tgen1, cX, cY)

%
% test 2: cX == 1 and cY > 1, no padding
%

cX = 1;
cY = randi(10, 1);

Tgen2 = randi(100, [cX, cY])
H2 = rooftopMoM.core.greensMatrix.utils.embedToeplitz4ToCirculant(Tgen2, cX, cY)

%
% test 3: cX > 1 and cY > 1, no padding
%

cX = randi(10, 1);
cY = randi(10, 1);

Tgen3 = randi(100, [cX, cY])
H3 = rooftopMoM.core.greensMatrix.utils.embedToeplitz4ToCirculant(Tgen3, cX, cY)

%
% test 4: cY == 1 and cX > 1, padding
%

cX = randi(10, 1);
cY = 1;

Tgen4 = randi(100, [cX, cY])
H4 = rooftopMoM.core.greensMatrix.utils.embedToeplitz4ToCirculant(Tgen4, cX + 1, cY)

%
% test 5: cX == 1 and cY > 1, padding
%

cX = 1;
cY = randi(10, 1);

Tgen5 = randi(100, [cX, cY])
H5 = rooftopMoM.core.greensMatrix.utils.embedToeplitz4ToCirculant(Tgen5, cX, cY + 1)

%
% test 6: cX > 1 and cY > 1, padding
%

cX = randi(10, 1);
cY = randi(10, 1);

Tgen6 = randi(100, [cX, cY])
H6 = rooftopMoM.core.greensMatrix.utils.embedToeplitz4ToCirculant(Tgen6, cX, cY + 1)

%
% test 7: cX = 0
%

% Results in an error.
H7 = rooftopMoM.core.greensMatrix.utils.embedToeplitz4ToCirculant(Tgen6, 0, cY);
