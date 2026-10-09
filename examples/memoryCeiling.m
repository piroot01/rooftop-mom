% memoryCeiling: The memory a run may take. A ceiling caps the data of rooftopMoM below the free memory of the machine;
%   check estimates a planned run before anything is built and stops it with what it needs if it does not fit. Every
%   object checks its own part again before it allocates.

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));

% The data may take 2 GB at most, and MATLAB 4 computational threads.
rooftopMoM.Resources.configure('memory', 2 * 2^30, 'threads', 4);
fprintf('Free memory %.1f GB, ceiling %.1f GB.\n\n', rooftopMoM.Resources.available() / 2^30, ...
    rooftopMoM.Resources.ceiling() / 2^30);

% A plate of 2048 x 2048 faces solved by GMRES does not fit; the error tells what it needs.
try
    rooftopMoM.Resources.check(rooftopMoM.Mesh(2048, 2048, 0.05, 0.05), 'method', 'gmres', 'restart', 50);
catch exception
    fprintf('%s\n\n', exception.message);
end

% A plate of 256 x 256 faces fits; the report gives the memory of every stage.
mesh = rooftopMoM.Mesh(256, 256, 0.05, 0.05);
report = rooftopMoM.Resources.check(mesh);
names = fieldnames(report);

for n = 1:numel(names)
    fprintf('%-22s %7.1f MB\n', names{n}, report.(names{n}) / 2^20);
end

problem = rooftopMoM.Problem(mesh, 2 * pi);
operator = rooftopMoM.Operator(problem);
preconditioner = rooftopMoM.Preconditioner(problem);
result = rooftopMoM.Solver().solve(operator, rooftopMoM.Excitation.planeWave(problem), ...
    'preconditioner', preconditioner);

fprintf('\nSolved in %d iterations; the live objects hold %.1f MB.\n', result.iterationCount(), ...
    rooftopMoM.Resources.held() / 2^20);

rooftopMoM.Resources.reset();
