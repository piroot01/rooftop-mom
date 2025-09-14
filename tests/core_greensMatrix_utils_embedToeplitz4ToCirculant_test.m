clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

%
% test 1: cY == 1 and cX > 1, no padding
%

cX = randi(10, 1);
cY = 1;

Tgen1 = randi(100, [cY, cX])
H1 = rooftopMoM.core.greensMatrix.utils.embedToeplitz4ToCirculant(Tgen1, cY, cX)

%
% test 2: cX == 1 and cY > 1, no padding
%

cX = 1;
cY = randi(10, 1);

Tgen2 = randi(100, [cY, cX])
H2 = rooftopMoM.core.greensMatrix.utils.embedToeplitz4ToCirculant(Tgen2, cY, cX)

%
% test 3: cX > 1 and cY > 1, no padding
%

cX = randi(10, 1);
cY = randi(10, 1);

Tgen3 = randi(100, [cY, cX])
H3 = rooftopMoM.core.greensMatrix.utils.embedToeplitz4ToCirculant(Tgen3, cY, cX)

%
% test 4: cY == 1 and cX > 1, padding
%

cX = randi(10, 1);
cY = 1;

Tgen4 = randi(100, [cY, cX])
H4 = rooftopMoM.core.greensMatrix.utils.embedToeplitz4ToCirculant(Tgen4, cY, cX + 1)

%
% test 5: cX == 1 and cY > 1, padding
%

cX = 1;
cY = randi(10, 1);

Tgen5 = randi(100, [cY, cX])
H5 = rooftopMoM.core.greensMatrix.utils.embedToeplitz4ToCirculant(Tgen5, cY + 1, cX)

%
% test 6: cX > 1 and cY > 1, padding
%

cX = randi(10, 1);
cY = randi(10, 1);

Tgen6 = randi(100, [cY, cX])
H6 = rooftopMoM.core.greensMatrix.utils.embedToeplitz4ToCirculant(Tgen6, cY + 1, cX)

%
% test 7: cX = 0
%

% Results in an error.
H7 = rooftopMoM.core.greensMatrix.utils.embedToeplitz4ToCirculant(Tgen6, cY, 0);
