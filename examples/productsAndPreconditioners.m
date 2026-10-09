% productsAndPreconditioners: The settings that trade time for memory, on a plate of 64 x 64 faces. The three
%   products give the same result at a different cost; the preconditioners differ in the iterations they save.

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));

% A finer quadrature than the default, and the memory of one assembly chunk set to 32 MB.
problem = rooftopMoM.Problem(rooftopMoM.Mesh(64, 64, 0.05, 0.05), 2 * pi, 'quadratureOrder', 6, ...
    'memoryBytes', 2^25);
input = rand(problem.rooftopCount(), 1);
products = {'fft', 'splitFFT', 'leanSplitFFT'};

fprintf('%d unknowns\n', problem.rooftopCount());

for n = 1:numel(products)
    operator = rooftopMoM.Operator(problem, 'product', products{n});
    fprintf('%-13s %6.2f ms per product, %5.1f MB stored, %5.1f MB at the peak\n', products{n}, ...
        1e3 * timeit(@() operator.apply(input)), operator.memory() / 2^20, operator.workingMemory() / 2^20);
end

operator = rooftopMoM.Operator(problem);
excitation = rooftopMoM.Excitation.planeWave(problem);
solver = rooftopMoM.Solver('tolerance', 1e-6);

fprintf('\n');
report_('no preconditioner', solver.solve(operator, excitation));

kernels = {'strang', 'tChan', 'rChan', 'huckle'};
builders = {'partwise', 'hybrid'};

for k = 1:numel(kernels)
    for b = 1:numel(builders)
        preconditioner = rooftopMoM.Preconditioner(problem, 'kernel', kernels{k}, 'builder', builders{b});
        report_([kernels{k}, ', ', builders{b}], solver.solve(operator, excitation, 'preconditioner', preconditioner));
    end
end

% The huckle kernel keeps the offsets within its bandwidth, by default the whole grid.
preconditioner = rooftopMoM.Preconditioner(problem, 'kernel', 'huckle', 'bandwidth', [16, 16]);
report_('huckle, bandwidth 16', solver.solve(operator, excitation, 'preconditioner', preconditioner));

% A modifier builds the preconditioner from a part of the matrix, here its imaginary part.
preconditioner = rooftopMoM.Preconditioner(problem, 'modifier', @imag);
report_('tChan, imaginary part', solver.solve(operator, excitation, 'preconditioner', preconditioner));

function report_(name, result)
    fprintf('%-22s %5d iterations %7.2f s\n', name, result.iterationCount(), result.time());
end
