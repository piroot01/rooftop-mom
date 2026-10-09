function bytes = estimateSolver(count, method, isPreconditioned, restart, maxIterationCount)
% estimateSolver: Returns the bytes of the vectors an iterative solve of count unknowns holds at its peak, as measured
%   in complex vectors of the system size:
%       'cg': the solution, the residual, the direction, its product and one temporary of an update, 5 vectors, and
%           the preconditioned residual with a preconditioner.
%       'gmres': the basis of a cycle, min(restart, count, maxIterationCount) + 1 vectors, half of it once more while
%           its fresh memory is first touched, 2 work vectors, and the preconditioned product with a preconditioner.
%
%   Inputs:
%       count: Number of unknowns.
%       method: 'cg' or 'gmres'.
%       isPreconditioned: Logical scalar.
%       restart: Iterations between the restarts of 'gmres'.
%       maxIterationCount: Maximum number of iterations.
%
%   Output:
%       bytes: Memory of the vectors.

    if strcmp(method, 'cg')
        vectorCount = 5 + isPreconditioned;
    else
        basisCount = min([restart, count, maxIterationCount]) + 1;
        vectorCount = 1.5 * basisCount + 2 + isPreconditioned;
    end

    bytes = 16 * count * vectorCount;
end
