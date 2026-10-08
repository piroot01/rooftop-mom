classdef (Abstract) IterativeMethod
% IterativeMethod: Base of the iterative solvers of A * x = b. A solve stops when the relative residual norm drops
%   below the tolerance or after the maximum number of iterations. The operator A (and a preconditioner) is a matrix,
%   a function handle or any object with the method apply.

    properties (Access = protected)
        tolerance_ double;
        maxIterationCount_ double;
    end

    methods
        function obj = IterativeMethod(tolerance, maxIterationCount)
            obj.tolerance_ = tolerance;
            obj.maxIterationCount_ = maxIterationCount;
        end
    end

    methods
        function tolerance = tolerance(obj)
            tolerance = obj.tolerance_;
        end

        function maxIterationCount = maxIterationCount(obj)
            maxIterationCount = obj.maxIterationCount_;
        end

        function output = solve(obj, operator, rightHandSide, initialGuess)
        % solve: Solves operator * x = rightHandSide.
        %
        %   Inputs:
        %       operator: Matrix, function handle or object with apply, the matrix A.
        %       rightHandSide: Nx1 vector b.
        %       initialGuess: Nx1 start vector (empty for zeros).
        %
        %   Output:
        %       output: rooftopMoM.solver.Output.
        %
        %   NOTE: The output holds the true residual norm(b - A * x) / norm(b) of the solution. It costs one more
        %   product, but the residual monitored by the method can drift away from it.

            start = tic;
            apply = rooftopMoM.solver.utils.convertToHandle(operator);
            logger = rooftopMoM.logging.Logger.instance();

            % The lines are written here for every method, so they name the method.
            method = regexprep(class(obj), '^.*\.', '');
            logger.debug('Solving %d unknowns with %s to the relative residual %.1e in at most %d iterations.', ...
                numel(rightHandSide), method, obj.tolerance_, obj.maxIterationCount_);

            if any(rightHandSide)
                [solution, converged, residualNorms] = obj.iterate_(apply, rightHandSide, initialGuess);
                residualNorm = norm(rightHandSide - apply(solution)) / norm(rightHandSide);
            else
                % The relative residual of a zero right-hand side is not defined; its solution is zero.
                solution = zeros(size(rightHandSide), 'like', rightHandSide);
                [converged, residualNorms, residualNorm] = deal(true, 0, 0);
            end

            iterationCount = numel(residualNorms) - 1;
            time = toc(start);
            output = rooftopMoM.solver.Output(solution, converged, iterationCount, residualNorms, residualNorm, time);

            if converged
                logger.info('%s converged in %d iterations in %.3f s, relative residual %.2e.', method, ...
                    iterationCount, time, residualNorm);
            else
                logger.warning('%s did not converge in %d iterations (%.3f s), relative residual %.2e.', method, ...
                    iterationCount, time, residualNorm);
            end
        end
    end

    methods (Abstract, Access = protected)
        % iterate_: Returns the solution, whether it converged, and the relative residual norms monitored from the
        %   initial guess to the last iteration. The right-hand side is not zero.
        [solution, converged, residualNorms] = iterate_(obj, apply, rightHandSide, initialGuess);
    end
end
