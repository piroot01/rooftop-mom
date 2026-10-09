% logToFile: Every stage of a run reports to one logger. By default only warnings reach the command window; a level
%   and a file make a run traceable, and the level trace adds a line per iteration and per chunk of the assembly. The
%   level debug also shows every memory check against the ceiling.

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));

logFile = fullfile(tempdir, 'rooftopMoM.log');

% The file is appended to, so the log of an earlier run goes first.
if isfile(logFile)
    delete(logFile);
end

rooftopMoM.Logging.configure('level', 'trace', 'file', logFile, 'console', false);
rooftopMoM.Resources.configure('memory', 2^28);

mesh = rooftopMoM.Mesh(16, 16, 0.05, 0.05);
rooftopMoM.Resources.check(mesh);
problem = rooftopMoM.Problem(mesh, 2 * pi);
rooftopMoM.Solver().solve(rooftopMoM.Operator(problem), rooftopMoM.Excitation.planeWave(problem), ...
    'preconditioner', rooftopMoM.Preconditioner(problem));

% Back to the defaults: no ceiling, warnings to the command window.
rooftopMoM.Resources.reset();
rooftopMoM.Logging.reset();
type(logFile);
