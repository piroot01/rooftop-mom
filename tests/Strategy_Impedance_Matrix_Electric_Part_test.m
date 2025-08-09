clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 3;
countY = 3;
sizeX = 1;
sizeY = 1;

mesh = rooftopMoM.mesh.utils.generateRectMesh(countX, countY, sizeX, sizeY);

settings = rooftopMoM.StrategySettings( ...
    rooftopMoM.StrategyOptions.IMPEDANCE_MATRIX_ELECTRIC_PART ...
);

settings.setOrder(5);
settings.setWavenumber(1);
settings.setElementCount(100000);

magStrat = rooftopMoM.Strategy(settings);

output = magStrat.generateWithMesh(mesh);

output.value()
