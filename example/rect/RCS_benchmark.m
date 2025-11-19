%% RCS benchmark
clc; clearvars;

% Include paths.
currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', '..', 'src');
addpath(srcDir);

utilsDir = fullfile(currentDir, 'utils');
addpath(utilsDir);

loggerDir = fullfile(currentDir, '..', '..', '3rd', 'matlab-logger', 'src');
addpath(loggerDir);

% Configure the logger.
logger.Logger.configure( ...
    'timeFormat', 'ss.SSS', ...
    'minLevel', logger.LogLevel.INFO, ...
    'sinks', { ...
        logger.sink.ConsoleSink(), ...
    } ...
);

output = solveAntennaCurrent();
output.solverOutput.postRunInfo().current
