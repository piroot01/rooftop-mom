function [lambda, V, info] = solveEig1Problem(eigCount, applicators, options)
% solve_eig1: Solves the imag(Z) * x = lambda * real(Z) * I using applicators.
%   Due to the bad convergence of solution I to equation real(Z) * I = V we tweak
%   the problem to an equivalent one imag(Z)^(-1) * real(Z) * I = (lambda)^(-1) * I.
%
%   Inputs:
%       eigCount: Number of requested eigenvalues.
%       applicators: Structs of the applicators associated with the problem.
%       options.innerTolerance (default 1e-6): Tolerance of the PCG solver.
%       options.innerMaxIterCount (default 1e3): Maximum number of iterations for the inner PCG solver.
%       options.outerTolerance (default 1e-6): Tolerance passed to eigs function.
%       options.outerMaxIterCount (default 1e3): Maximum number of iterations for the eigs function.
%       options.pcgWarn (default true): Shows warnings if the PCG fails.
%       options.whichEigs (default 'largestabs'): Eigenvalues which to search.
%
%   Outputs:
%       lambda: 1 by k vector of eigenvalues.
%       V: n by k matrix of eigenvector.
%       info: Structure with useful informations about the solution run.

    if nargin < 5
        options = struct();
    end

    innerTolerance = getOption(options, 'innerTolerance', 1e-6);
    innerMaxIterCount = getOption(options, 'innerMaxIterCount', 1e3);
    outerTolerance = getOption(options, 'outerTolerance', 1e-6);
    outerMaxIterCount = getOption(options, 'outerMaxIterCount', 1e3);
    pcgWarn = getOption(options, 'pcgWarn', true);
    whichEigs = getOption(options, 'whichEigs', 'largestabs');

    info.innerCallsCount = 0;
    info.pcgFailCount = 0;
    info.pcgTotalIterCount = 0;

    % Inner solver that solves imag(Z) * I = V for I, we use PCG
    function output = applyInvImagZ(solver, input)
        solverOutput = solver.solve( ...
            'initCurrent', zeros(applicators.imagGreensMatrix.inputSize(), 1), ...
            'excitation', input, ...
            'tolerance', innerTolerance, ...
            'maxIterCount', innerMaxIterCount ...
        );

        output = complex(solverOutput.postRunInfo().current);

        info.pcgTotalIterCount = info.pcgTotalIterCount + solverOutput.count();

        if ~solverOutput.postRunInfo().converged
            if pcgWarn
                WARN('Inner PCG solver does not converged.');
            end

            info.pcgFailCount = info.pcgFailCount + 1;
        end
    end

    % Define the composed operator of the equivalent eigenvalue problem
    function output = applyInvImagZRealZ(solver, input)
        [nrows, ncols] = size(input);

        output = zeros(nrows, ncols, 'like', input);

        for col = 1:ncols
            t = applicators.realGreensMatrix.apply(input(:, col));
            output(:, col) = applyInvImagZ(solver, t);
        end
    end

    % Initialize the solver for the inner part.
    solver = rooftopMoM.core.solver.PreconditionedConjugateGradient( ...
        applicators.imagGreensMatrix, applicators.imagPreconditioner);

    % Make a function handler for the eigs function.
    eigsHandle = @(x) applyInvImagZRealZ(solver, x);

    [V, D] = eigs(eigsHandle, ...
        applicators.realGreensMatrix.inputSize(), ...
        eigCount, whichEigs, ...
        struct('tol', outerTolerance, 'maxit', outerMaxIterCount));

    % And finally invert.
    lambda = 1 ./ diag(D);
end
