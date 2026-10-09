function measureFFTPlanner(count, planners)
% measureFFTPlanner: Prints the time of the three products and the circulant preconditioner of a plate of count x count
%   faces of a tenth of a wavelength for each planner method of FFTW. A slower planner tries more algorithms once per
%   transform size and may give faster transforms; the set-up column shows that cost, the other columns the medians
%   of the products after it.
%       measureFFTPlanner(2048, {'estimate', 'measure', 'patient'})
%   Every method starts without wisdom, so it plans for itself. For a strict comparison, run each method in its own
%   MATLAB session. A method that pays off can keep its plans across sessions:
%       wisdom = fftw('dwisdom');     % after the products ran, saved with save()
%       fftw('dwisdom', wisdom);      % in the next session, before the products
%
%   Inputs:
%       count: Number of faces along each side (default 1024).
%       planners: Cell array of the methods (default {'estimate', 'measure', 'patient'}).

    addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));
    rooftopMoM.Resources.reset();

    if nargin < 1
        count = 1024;
    end

    if nargin < 2
        planners = {'estimate', 'measure', 'patient'};
    end

    problem = rooftopMoM.Problem(rooftopMoM.Mesh(count, count, 0.1, 0.1), 2 * pi);
    problem.generators();
    input = complex(randn(problem.rooftopCount(), 1), randn(problem.rooftopCount(), 1));
    names = {'fft', 'splitFFT', 'leanSplitFFT', 'preconditioner'};

    fprintf('%d x %d faces, %d unknowns, %d threads, MATLAB %s\n\n%10s', count, count, problem.rooftopCount(), ...
        maxNumCompThreads, version('-release'), 'planner');
    fprintf('%25s', names{:});
    fprintf('\n');
    previous = fftw('planner');

    for k = 1:numel(planners)
        fftw('dwisdom', []);
        fftw('planner', planners{k});
        fprintf('%10s', planners{k});

        for n = 1:numel(names)
            % The set-up builds the spectra and plans the transforms of the first product.
            start = tic;

            if n < numel(names)
                product = rooftopMoM.Operator(problem, 'product', names{n});
            else
                product = rooftopMoM.Preconditioner(problem);
            end

            product.apply(input);
            setup = toc(start);
            times = zeros(1, 5);

            for m = 1:numel(times)
                start = tic;
                product.apply(input);
                times(m) = toc(start);
            end

            fprintf('%10.2f s + %8.3f s', setup, median(times));
            clear product
        end

        fprintf('\n');
    end

    fftw('planner', previous);
end
