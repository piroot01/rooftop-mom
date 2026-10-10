% productsAndPreconditioners: The settings that trade time for memory, on a plate of 64 x 64 faces. The three
%   products give the same result at a different cost; the preconditioners differ in the iterations they save.

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));

% The comparison may take 512 MB, and one assembly chunk 32 MB. The quadrature orders [12, 6] of the near and the far
% cell pairs assemble faster than the default [16, 8] and still give the entries to about 1e-12.
rooftopMoM.Resources.configure('memory', 2^29, 'chunkSize', 2^25);
mesh = rooftopMoM.Mesh(64, 64, 0.05, 0.05);
problem = rooftopMoM.Problem(mesh, 2 * pi, 'quadratureOrder', [12, 6]);
input = rand(problem.rooftopCount(), 1);
products = {'fft', 'splitFFT', 'leanSplitFFT'};

fprintf('%d unknowns\n', problem.rooftopCount());

for n = 1:numel(products)
    % check estimates the run with each product before the product is built.
    report = rooftopMoM.Resources.check(mesh, 'product', products{n});
    operator = rooftopMoM.Operator(problem, 'product', products{n});
    fprintf('%-13s %6.2f ms per product, %5.1f MB stored (%5.1f MB estimated), %5.1f MB at the peak\n', ...
        products{n}, 1e3 * timeit(@() operator.apply(input)), operator.memory() / 2^20, report.product / 2^20, ...
        operator.workingMemory() / 2^20);
end

operator = rooftopMoM.Operator(problem);
excitation = rooftopMoM.Excitation.planeWave(problem);
solver = rooftopMoM.Solver('tolerance', 1e-6);

fprintf('\n');
report_('no preconditioner', solver.solve(operator, excitation));

kernels = {'strang', 'tChan', 'rChan', 'huckle'};
builders = {'partwise', 'hybrid'};

for k = 1:numel(kernels)
    for n = 1:numel(builders)
        preconditioner = rooftopMoM.Preconditioner(problem, 'kernel', kernels{k}, 'builder', builders{n});
        report_([kernels{k}, ', ', builders{n}], solver.solve(operator, excitation, 'preconditioner', preconditioner));
    end
end

% The blockwise builder approximates the four blocks of the matrix in space; it suits the R. Chan kernel.
preconditioner = rooftopMoM.Preconditioner(problem, 'kernel', 'rChan', 'builder', 'blockwise');
report_('rChan, blockwise', solver.solve(operator, excitation, 'preconditioner', preconditioner));

% The superoptimal circulant of Tyrtyshnikov, a kernel like the others.
preconditioner = rooftopMoM.Preconditioner(problem, 'kernel', 'superoptimal');
report_('superoptimal, partwise', solver.solve(operator, excitation, 'preconditioner', preconditioner));

% The huckle kernel keeps the offsets within its bandwidth, by default the whole grid.
preconditioner = rooftopMoM.Preconditioner(problem, 'kernel', 'huckle', 'bandwidth', [16, 16]);
report_('huckle, bandwidth 16', solver.solve(operator, excitation, 'preconditioner', preconditioner));

% A modifier builds the preconditioner from a part of the matrix, here its imaginary part.
preconditioner = rooftopMoM.Preconditioner(problem, 'modifier', @imag);
report_('tChan, imaginary part', solver.solve(operator, excitation, 'preconditioner', preconditioner));

% The generators, the product and the last preconditioner are still alive.
fprintf('\nThe live objects hold %.1f MB.\n', rooftopMoM.Resources.held() / 2^20);
rooftopMoM.Resources.reset();

function report_(name, result)
    fprintf('%-22s %5d iterations %7.2f s\n', name, result.iterationCount(), result.time());
end
