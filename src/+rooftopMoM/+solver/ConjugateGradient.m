classdef ConjugateGradient < rooftopMoM.solver.IterativeMethod
% ConjugateGradient: Conjugate gradient method for complex symmetric systems (COCG): the transpose replaces the
%   conjugate transpose in the inner products. It monitors the relative norm of the residual, norm(r) / norm(b).

    methods
        function obj = ConjugateGradient(tolerance, maxIterationCount)
            obj@rooftopMoM.solver.IterativeMethod(tolerance, maxIterationCount);
        end
    end

    methods (Access = protected)
        function [solution, converged, residualNorms] = iterate_(obj, apply, rightHandSide, initialGuess)
        % iterate_: Runs the COCG iteration, see rooftopMoM.solver.IterativeMethod.

            logger = rooftopMoM.logging.Logger.instance();
            isTrace = logger.isEnabled(rooftopMoM.logging.LogLevel.TRACE);
            rightHandSideNorm = norm(rightHandSide);

            if isempty(initialGuess)
                solution = zeros(size(rightHandSide), 'like', rightHandSide);
                residual = rightHandSide;
            else
                solution = initialGuess;
                residual = rightHandSide - apply(solution);
            end

            residualNorms = zeros(obj.maxIterationCount_ + 1, 1);
            residualNorms(1) = norm(residual) / rightHandSideNorm;
            converged = residualNorms(1) < obj.tolerance_;
            iteration = 0;

            while ~converged && iteration < obj.maxIterationCount_
                iteration = iteration + 1;
                preconditioned = obj.precondition_(residual);
                rhoNext = residual.' * preconditioned;

                if iteration == 1
                    direction = preconditioned;
                else
                    direction = preconditioned + (rhoNext / rho) * direction;
                end

                rho = rhoNext;
                appliedDirection = apply(direction);
                alpha = rho / (direction.' * appliedDirection);
                solution = solution + alpha * direction;
                residual = residual - alpha * appliedDirection;
                residualNorms(iteration + 1) = norm(residual) / rightHandSideNorm;
                converged = residualNorms(iteration + 1) < obj.tolerance_;

                if isTrace
                    logger.trace('Iteration %d: relative residual %.3e.', iteration, residualNorms(iteration + 1));
                end

                % A breakdown (direction.' * A * direction = 0) leaves Inf or NaN, so further iterations are wasted.
                if ~isfinite(residualNorms(iteration + 1))
                    break;
                end
            end

            residualNorms = residualNorms(1:(iteration + 1));
        end

        function preconditioned = precondition_(~, residual)
        % precondition_: Returns the preconditioned residual; the plain method has no preconditioner.

            preconditioned = residual;
        end
    end
end
