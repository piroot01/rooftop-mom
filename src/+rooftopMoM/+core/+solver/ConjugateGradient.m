classdef ConjugateGradient < handle
% ConjugateGradient: MATLAB implementation of a (matrix-free) Conjugate Gradient solver

    properties (Access = protected)
        matrixApplicator_;
    end

    methods
        function obj = ConjugateGradient(matrixApplicator)
            obj.matrixApplicator_ = matrixApplicator;
        end
    end

    methods
        function output = solve(obj, varargin)
        % solve: High-level entry point for the solver.
        %
        %   Optional name-value inputs:
        %       'initCurrent': Initial guess vector (complex, column).
        %       'excitation': Right-hand side vector (required).
        %       'tolerance': Relative residual tolerance (default 1e-6).
        %       'maxIterCount': Maximum number of iterations (default 1000).

            TRACE("ConjugateGradient.solve: Initialized solving procedure.");

            p = inputParser();

            addParameter(p, 'initCurrent', []);
            addParameter(p, 'excitation', []);
            addParameter(p, 'tolerance', 1e-6);
            addParameter(p, 'maxIterCount', 1000);

            parse(p, varargin{:});

            if isempty(p.Results.excitation)
                error('The excitation is empty!');
            end

            DEBUG("ConjugateGradient.solve: Parsed inputs: rooftop count = %d, tolerance = %g, maxIterCount = %d", ...
                numel(p.Results.initCurrent), p.Results.tolerance, p.Results.maxIterCount);

            output = obj.solveImpl_( ...
                p.Results.initCurrent, ...
                p.Results.excitation, ...
                p.Results.tolerance, ...
                p.Results.maxIterCount ...
            );
        end
    end

    methods (Access = protected)
        function output = solveImpl_(obj, initCurrent, excitation, tolerance, maxIterCount)
        % solveImpl_: Core Conjugate-Gradient loop.
        %
        %   Inputs:
        %       initCurrent: Initial guess column vector.
        %       excitation: Right-hand side column vector.
        %       tolerance: Relative tolerance.
        %       maxIterCount: Maximum iterations.
        %
        %   Returns:
        %       output: Solution in form of the ConjugateGradientOutput.

            TRACE("ConjugateGradient.solveImpl_: Starting CG solver core.");

            if isa(obj.matrixApplicator_, 'function_handle')
                applyHandle = obj.matrixApplicator_;
            else
                applyHandle = @(x) obj.matrixApplicator_.apply(x);
            end

            output = rooftopMoM.core.solver.ConjugateGradientOutput();

            current = initCurrent;
            residue = excitation - applyHandle(current);
            direction = residue;
            rTr = residue.' * residue;
            excitationNorm = norm(excitation);
            converged = false;

            DEBUG("ConjugateGradient.solveImpl_: Initial residual norm = %g, excitationNorm = %g.", sqrt(rTr) / excitationNorm, excitationNorm);

            t0 = tic;

            for iteration = 1:maxIterCount
                t1 = tic;

                gDirection = applyHandle(direction);
                alpha = rTr / (direction.' * gDirection);
                current = current + alpha * direction;
                residue = residue - alpha * gDirection;
                rTrNew = residue.' * residue;
                residueNorm = sqrt(rTrNew) / excitationNorm;

                if residueNorm < tolerance
                    DEBUG("ConjugateGradient.solveImpl_: Converged in %d iterations (residue norm = %g).", iteration, residueNorm);
                    converged = true;
                    break;
                end

                beta = rTrNew / rTr;
                direction = residue + beta * direction;
                rTr = rTrNew;

                % Insert iteration info.
                output.insertIterationInfo( ...
                    'residueNorm', residueNorm ...
                );

                iterationTime = toc(t1);
                TRACE("ConjugateGradient.solveImpl_: Iteration %d: residual norm = %g, iteration time = %g s.", iteration, residueNorm, iterationTime);
            end

            solvingTime = toc(t0);

            if ~converged
                DEBUG("ConjugateGradient.solveImpl_: Did not converge after %d iterations (final residue norm = %g).", maxIterCount, residueNorm);
            end

            DEBUG("ConjugateGradient.solveImpl_: Total time = %g s.", solvingTime);

            output.setPostRunInfo( ...
                struct( ...
                    'converged', converged, ...
                    'current', current, ...
                    'solvingTime', solvingTime, ...
                    'residueNorm', residueNorm ...
                ) ...
            );

            TRACE("ConjugateGradient.solveImpl_: Exiting CG solver core.");
        end
    end
end
