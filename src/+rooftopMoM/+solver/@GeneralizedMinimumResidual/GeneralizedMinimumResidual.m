classdef GeneralizedMinimumResidual < rooftopMoM.solver.IterativeMethod
% GeneralizedMinimumResidual: Restarted GMRES for general systems, with an optional preconditioner P applied from the
%   left, P * A * x = P * b. It monitors the residual of the preconditioned system relative to P * b. A cycle builds an
%   orthonormal basis of up to restart + 1 vectors by two passes of Gram-Schmidt and solves its small least squares
%   problem by Givens rotations.

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
        [solution, converged, residualNorms] = iterate_(obj, apply, rightHandSide, initialGuess);
    end

    methods (Access = private)
        output = precondition_(obj, input);
    end
end
