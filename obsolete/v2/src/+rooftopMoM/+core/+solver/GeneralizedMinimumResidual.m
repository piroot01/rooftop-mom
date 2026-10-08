classdef GeneralizedMinimumResidual < rooftopMoM.core.solver.GreensMatrixApplicator
% GeneralizedMinimumResidual: Wrapper using MATLAB's built-in gmres with

    methods
        function obj = GeneralizedMinimumResidual(context, greensMatrix)
            obj@rooftopMoM.core.solver.GreensMatrixApplicator(context, greensMatrix);
        end
    end

    methods
        function output = solve(obj, varargin)
        % solve: High-level entry point for the solver.
        %
        %   Optional name-value inputs:
        %       'initCurrent': initial guess (complex column) [default zeros]
        %       'excitation': RHS (required)
        %       'tolerance': relative tol [default 1e-6]
        %       'maxIterCount': total maximum inner iterations (default 1000)
        %       'restart': GMRES restart dimension (default 50)

            TRACE("GeneralizedMinimumResidual.solve: initializing.");

            count = obj.context_.rooftops().count();

            p = inputParser();
            addParameter(p, 'initCurrent', complex(zeros(count, 1)));
            addParameter(p, 'excitation', []);
            addParameter(p, 'tolerance', 1e-6);
            addParameter(p, 'maxIterCount', 1000);
            addParameter(p, 'restart', 50);
            parse(p, varargin{:});

            if isempty(p.Results.excitation)
                error('The excitation is empty!');
            end

            DEBUG("GeneralizedMinimumResidual.solve: Parsed inputs: rooftop count = %d, tolerance = %g, maxIterCount = %d", ...
                count, p.Results.tolerance, p.Results.maxIterCount);

            output = obj.solveImpl_( ...
                p.Results.initCurrent, ...
                p.Results.excitation, ...
                p.Results.tolerance, ...
                p.Results.maxIterCount, ...
                p.Results.restart ...
            );
        end
    end

    methods (Access = protected)
        function output = solveImpl_(obj, initCurrent, excitation, tolerance, maxIterCount, restart)
            TRACE("GeneralizedMinimumResidual.solveImpl_: Starting GMRES solver core.");

            output = rooftopMoM.core.solver.GeneralizedMinimumResidualOutput();

            % Create the matrix handle, in our case the custom applyGreensMatrix_ function represents the matrix multiplication.
            matrixHandle = @(x) obj.applyGreensMatrix_(x);

            % Map user's maxIterCount (total inner iters) to gmres outer iters.
            % MATLAB gmres: gmres(Afun,b,restart,tol,maxit,...) where maxit = number of outer iterations.
            outerMaxIt = max(1, ceil(double(maxIterCount) / double(restart)));

            % Run gmres and request residual history (resvec).
            t0 = tic;

            try
                [x, flag, relres, iter, resvec] = gmres(matrixHandle, excitation, restart, tolerance, outerMaxIt, [], [], initCurrent);
            catch ME
                % If gmres errors, still record and rethrow.
                solvingTime = toc(t0); %#ok

                DEBUG("GeneralizedMinimumResidual.solveImpl_: GMRES failed: %s", ME.message);

                rethrow(ME);
            end

            solvingTime = toc(t0);

            DEBUG("GeneralizedMinimumResidual.solveImpl_: Total time = %g s.", solvingTime);

            % First element is initial residual norm (||b - A*x0||).
            excitationNorm = norm(excitation);
            if excitationNorm == 0
                excitationNorm = 1;
            end

            % Insert iteration info for each inner iteration.
            totalInnerPerformed = numel(resvec) - 1;
            for k = 2:numel(resvec)
                % Relative residual after that inner interation.
                rel = resvec(k) / excitationNorm;
                output.insertIterationInfo('residueNorm', rel);
            end

            % Post run info (consistent with ConjugateGradientOutput)
            converged = (flag == 0);

            % Build final info struct
            output.setPostRunInfo( ...
                struct( ...
                    'converged', converged, ...
                    'current', x, ...
                    'solvingTime', solvingTime, ...
                    'residueNorm', (resvec(end) / excitationNorm), ...
                    'totalInnerIterations', totalInnerPerformed, ...
                    'gmresFlag', flag, ...
                    'gmresRelRes', relres, ...
                    'gmresIter', iter ...
                ) ...
            );

            TRACE("GeneralizedMinimumResidual.solveImpl_: GMRES finished (flag=%d, relres=%g).", flag, relres);
        end
    end
end
