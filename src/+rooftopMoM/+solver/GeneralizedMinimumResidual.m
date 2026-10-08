classdef GeneralizedMinimumResidual < rooftopMoM.solver.IterativeMethod
% GeneralizedMinimumResidual: Restarted GMRES of MATLAB (gmres) for general systems, with an optional preconditioner.
%   The iteration limit counts single iterations and is rounded up to whole restart cycles. The monitored residual is
%   that of the preconditioned system.

    properties (Access = private)
        restart_ double;

        % Function handle of the preconditioner, empty for none.
        preconditioner_;
    end

    methods
        function obj = GeneralizedMinimumResidual(tolerance, maxIterationCount, restart, preconditioner)
            obj@rooftopMoM.solver.IterativeMethod(tolerance, maxIterationCount);
            obj.restart_ = restart;

            if ~isempty(preconditioner)
                obj.preconditioner_ = rooftopMoM.solver.utils.convertToHandle(preconditioner);
            end
        end
    end

    methods
        function restart = restart(obj)
            restart = obj.restart_;
        end
    end

    methods (Access = protected)
        function [solution, converged, residualNorms] = iterate_(obj, apply, rightHandSide, initialGuess)
        % iterate_: Runs gmres, see rooftopMoM.solver.IterativeMethod.

            count = numel(rightHandSide);

            % Without restarts (an empty restart) the limit of gmres counts single iterations, at most count of them.
            if obj.restart_ < count
                restart = obj.restart_;
                limit = ceil(obj.maxIterationCount_ / restart);
            else
                restart = [];
                limit = min(obj.maxIterationCount_, count);
            end

            % The residual vector holds the initial residual and one value per iteration.
            [solution, flag, ~, ~, residualVector] = gmres(apply, rightHandSide, restart, obj.tolerance_, limit, ...
                obj.preconditioner_, [], initialGuess);
            converged = flag == 0;

            % gmres monitors M \ (b - A * x), relative to M \ b.
            if isempty(obj.preconditioner_)
                residualNorms = residualVector / norm(rightHandSide);
            else
                residualNorms = residualVector / norm(obj.preconditioner_(rightHandSide));
            end
        end
    end
end
