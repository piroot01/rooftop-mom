classdef PreconditionedConjugateGradient < rooftopMoM.solver.ConjugateGradient
% PreconditionedConjugateGradient: Conjugate gradient method for complex symmetric systems (COCG) with a preconditioner,
%   an approximate inverse of A applied to the residual in every iteration.

    properties (Access = private)
        % Function handle of the preconditioner.
        preconditioner_;
    end

    methods
        function obj = PreconditionedConjugateGradient(preconditioner, tolerance, maxIterationCount)
            obj@rooftopMoM.solver.ConjugateGradient(tolerance, maxIterationCount);
            obj.preconditioner_ = rooftopMoM.solver.utils.convertToHandle(preconditioner);
        end
    end

    methods (Access = protected)
        function preconditioned = precondition_(obj, residual)
        % precondition_: Returns the preconditioner applied to the residual.

            preconditioned = obj.preconditioner_(residual);
        end
    end
end
