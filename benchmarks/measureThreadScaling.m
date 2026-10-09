function measureThreadScaling(count, threads)
% measureThreadScaling: Prints the wall time and the busy cores of the assembly, the three products and the circulant
%   preconditioner of a plate of count x count faces of a tenth of a wavelength, for each number of threads.
%       measureThreadScaling(2048, [1, 2, 4, 8, 16])
%   The times are medians after one warm-up run; the busy cores are the CPU time of MATLAB over the wall time. A run
%   that keeps many cores busy without getting faster waits for the memory, not for the cores.
%
%   Inputs:
%       count: Number of faces along each side (default 1024).
%       threads: Numbers of computational threads (default 1, 2, 4, ... and all cores).

    addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));
    rooftopMoM.Resources.reset();

    if nargin < 1
        count = 1024;
    end

    if nargin < 2
        threads = unique([2.^(0:floor(log2(feature('numcores')))), feature('numcores')]);
    end

    mesh = rooftopMoM.Mesh(count, count, 0.1, 0.1);
    names = {'assembly', 'fft', 'splitFFT', 'leanSplitFFT', 'preconditioner'};
    fprintf('%s\n%d x %d faces, %d unknowns\n\n%8s', describeMachine_(), count, count, 2 * count * (count - 1), ...
        'threads');
    fprintf('%28s', names{:});
    fprintf('\n');

    for n = threads
        rooftopMoM.Resources.configure('threads', n);
        problem = rooftopMoM.Problem(mesh, 2 * pi);
        input = complex(randn(problem.rooftopCount(), 1), randn(problem.rooftopCount(), 1));
        [wall, busy] = deal(zeros(1, numel(names)));

        % A new problem assembles its generators again.
        [wall(1), busy(1)] = measure_(@() rooftopMoM.Problem(mesh, 2 * pi).generators(), 3);

        for k = 2:4
            operator = rooftopMoM.Operator(problem, 'product', names{k});
            [wall(k), busy(k)] = measure_(@() operator.apply(input), 5);
            clear operator
        end

        preconditioner = rooftopMoM.Preconditioner(problem);
        [wall(5), busy(5)] = measure_(@() preconditioner.apply(input), 5);
        clear preconditioner problem

        fprintf('%8d', n);
        fprintf('%14.3f s (%4.1f busy)', [wall; busy]);
        fprintf('\n');
    end

    rooftopMoM.Resources.reset();
end

function [wall, busy] = measure_(operation, repeats)
% measure_: Returns the median wall time of the operation after one warm-up run, and the mean number of busy cores.

    operation();
    times = zeros(1, repeats);
    cpu = cputime;
    total = tic;

    for k = 1:repeats
        start = tic;
        operation();
        times(k) = toc(start);
    end

    busy = (cputime - cpu) / toc(total);
    wall = median(times);
end

function text = describeMachine_()
% describeMachine_: Returns the processor, the cores and the MATLAB release of this machine.

    if ispc
        processor = getenv('PROCESSOR_IDENTIFIER');
    elseif ismac
        [~, processor] = system('sysctl -n machdep.cpu.brand_string');
    else
        processor = regexp(fileread('/proc/cpuinfo'), 'model name\s*:\s*([^\n]*)', 'tokens', 'once');
        processor = processor{1};
    end

    text = sprintf('%s, %d cores, MATLAB %s', strtrim(processor), feature('numcores'), version('-release'));
end
