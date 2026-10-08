classdef Output
% Output: Result of an iterative solve of A * x = b.

    properties (Access = private)
        solution_;
        converged_ logical;
        iterationCount_ double;
        residualNorms_ double;
        residualNorm_ double;
        time_ double;
    end

    methods
        function obj = Output(solution, converged, iterationCount, residualNorms, residualNorm, time)
            obj.solution_ = solution;
            obj.converged_ = converged;
            obj.iterationCount_ = iterationCount;
            obj.residualNorms_ = residualNorms;
            obj.residualNorm_ = residualNorm;
            obj.time_ = time;
        end
    end

    methods
        function solution = solution(obj)
            solution = obj.solution_;
        end

        function converged = converged(obj)
            converged = obj.converged_;
        end

        function iterationCount = iterationCount(obj)
            iterationCount = obj.iterationCount_;
        end

        function residualNorms = residualNorms(obj)
        % residualNorms: Returns the relative residual norms monitored by the method, from the initial guess to the
        %   last iteration (iterationCount + 1 values).

            residualNorms = obj.residualNorms_;
        end

        function residualNorm = residualNorm(obj)
        % residualNorm: Returns the relative norm of the true residual, norm(b - A * x) / norm(b).

            residualNorm = obj.residualNorm_;
        end

        function time = time(obj)
        % time: Returns the time of the solve in seconds.

            time = obj.time_;
        end
    end
end
