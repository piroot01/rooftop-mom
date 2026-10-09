classdef (TestTags = {'Slow'}) EstimateSolverTest < unit.solver.Base
% EstimateSolverTest: The estimate of rooftopMoM.resources.estimateSolver against the peak resident memory of real
%   solves. The vectors have 4 Mi complex unknowns, 64 MiB each, so every one of them is mapped and unmapped on its
%   own and the peak counts the live vectors only.

    properties (TestParameter)
        setting = struct( ...
            'cg', struct('method', 'cg', 'isPreconditioned', false, 'restart', Inf), ...
            'preconditionedCg', struct('method', 'cg', 'isPreconditioned', true, 'restart', Inf), ...
            'gmres', struct('method', 'gmres', 'isPreconditioned', false, 'restart', 10), ...
            'preconditionedGmres', struct('method', 'gmres', 'isPreconditioned', true, 'restart', 10));
    end

    methods (Test)
        function testCoversThePeakOfASolve(testCase, setting)
            % Linux keeps the peak resident memory of the process and resets it on request.
            testCase.assumeTrue(isfile('/proc/self/clear_refs'), 'The peak memory is read from Linux /proc only.');

            count = 2^22;
            diagonal = complex(linspace(1, 50, count).', 1);
            rightHandSide = complex(ones(count, 1), 0);
            preconditioner = [];

            if setting.isPreconditioned
                preconditioner = @(residual) residual ./ sqrt(diagonal);
            end

            method = createMethod_(setting, preconditioner);

            % The first large solve of a session costs more once; a solve of the same size before settles it.
            method.solve(@(input) diagonal .* input, rightHandSide, []);
            baseline = testCase.readProcessMemory('VmRSS');
            testCase.resetPeakMemory();
            method.solve(@(input) diagonal .* input, rightHandSide, []);
            peak = testCase.readProcessMemory('VmHWM') - baseline;

            estimate = rooftopMoM.resources.estimateSolver(count, setting.method, setting.isPreconditioned, ...
                setting.restart, 30);

            % The estimate of GMRES holds half its basis once more for the first touch of fresh memory, which not every
            % solve shows; without it the peak is about 0.72 of the estimate.
            testCase.verifyLessThanOrEqual(peak, estimate);
            testCase.verifyGreaterThanOrEqual(peak, 0.7 * estimate);
        end
    end
end

function method = createMethod_(setting, preconditioner)
    if strcmp(setting.method, 'cg')
        if isempty(preconditioner)
            method = rooftopMoM.solver.ConjugateGradient(1e-14, 30);
        else
            method = rooftopMoM.solver.PreconditionedConjugateGradient(preconditioner, 1e-14, 30);
        end
    else
        method = rooftopMoM.solver.GeneralizedMinimumResidual(1e-14, 30, setting.restart, preconditioner);
    end
end
