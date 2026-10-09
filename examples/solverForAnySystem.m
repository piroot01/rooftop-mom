% solverForAnySystem: The solver needs only its own module, the memory budget and the logging of the package; it takes
%   a matrix, a function handle or any object with the method apply. Here it solves the Helmholtz equation with a
%   little loss on the unit square, whose finite-difference matrix is complex symmetric like the Green's matrix.

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));

% Before a solve the solver checks its vectors against the ceiling; one core stays free for the rest of the machine.
rooftopMoM.Resources.configure('memory', 2^27, 'threads', max(1, maxNumCompThreads - 1));

% The five-point Laplacian on 100 x 100 inner points, minus the squared wavenumber.
count = 100;
wavenumber = 20 - 1i;
second = (count + 1)^2 * spdiags(ones(count, 1) * [-1, 2, -1], -1:1, count, count);
laplacian = kron(speye(count), second) + kron(second, speye(count));
matrix = laplacian - wavenumber^2 * speye(count^2);

% A point source in the middle.
rightHandSide = zeros(count^2, 1);
rightHandSide(sub2ind([count, count], count / 2, count / 2)) = (count + 1)^2;

% COCG on the matrix itself.
plain = rooftopMoM.Solver('maxIterationCount', 5000).solve(matrix, rightHandSide);

% GMRES on a function handle, preconditioned by a shifted Laplacian factorized once.
shifted = decomposition(laplacian - (1 - 0.5i) * real(wavenumber)^2 * speye(count^2));
preconditioned = rooftopMoM.Solver('method', 'gmres', 'restart', 30).solve(@(input) matrix * input, ...
    rightHandSide, 'preconditioner', @(residual) shifted \ residual);

% A basis of 1001 vectors does not fit below the ceiling, so this solve stops before it starts.
try
    rooftopMoM.Solver('method', 'gmres', 'restart', 1000).solve(matrix, rightHandSide);
catch exception
    disp(exception.message);
end

rooftopMoM.Resources.reset();

semilogy(0:plain.iterationCount(), plain.residualNorms(), ...
    0:preconditioned.iterationCount(), preconditioned.residualNorms());
xlabel('iteration');
ylabel('relative residual norm');
legend('COCG', 'GMRES with a shifted Laplacian');
grid on;
