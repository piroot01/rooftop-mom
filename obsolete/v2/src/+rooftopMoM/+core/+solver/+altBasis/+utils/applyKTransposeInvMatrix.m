function I = applyKTransposeInvMatrix(Q, treePatches)
% applyKTransposeInvMatrix: Solves a system of linear equations written in the matrix form K^T * I = Q.
%   The equation can be solved by applying the inverse of K^T to the right-hand side: I = (K^T)^(-1) * Q.
%   By leveraging the known structure of the matrix, we can accelerate the application of (K^T)^(-1)
%   substantially, achieving linear time complexity and essentially constant memory.
%
%   Inputs:
%       Q: Nx1 vector defining the right-hand side of the equation.
%       treePatches: Nx1 or 1xN vector defining how the tree spans the patches, and in what order.
%
%   Output:
%       I: Nx1 vector that solves the equation for the given Q and the K^T matrix defined by the tree.
%
%   Notes:
%       The K^T matrix is defined as the transpose of the Gram matrix described in the applyKInvMatrix function.
%
%   TODO: Currently, the function does not support rooftops with edges of non-unit length. This is
%       because the Haar function would then depend on the edge length of the rooftop inside the
%       inner product. Specifically, div(R^{(T)_{n}}(r)) = 1 / l_n for r on the positive side of the 
%       rooftop and -1 / l_n for r on the negative side.

    N = numel(Q);
    I = zeros(N, 1, like = Q);

    if N == 0
        return;
    end

    if N == 1
        % If the matrix is 1x1, I(1) simply equals Q(1).
        I(1) = Q(1);
        return;
    end

    % Pick the relevant pivot indices from the treePatches.
    pivots = treePatches(2:(end - 1));

    % Find the first occurrence where a pivot index exceeds N. This essentially marks the column
    % of the matrix where the reverse forward substitution transitions into a reverse backward substitution.
    idx = find(pivots > N, 1);

    if isempty(idx)
        % No missing row; forward substitution covers the whole array.
        idx = N;
    end

    % Index marking the first pivot of the backward substitution.
    dirChangeIdx = idx + 1;

    % Reverse Forward Substitution:
    if idx > 0
        F = flipud(cumsum(flipud(Q(1:idx))));

        I(1) = F(1);

        % The following part is only relevant if the index of the substitution change is greater than 1.
        if idx > 1
            forwardPivots = pivots(1:(idx - 1));
            I(forwardPivots) = I(forwardPivots) + F(2:end);
        end
    end

    % Forward cumulative sum applied backward:
    if dirChangeIdx <= N
        B = cumsum(Q(dirChangeIdx:N));

        % Determine the backward pivot indices.
        backwardPivots = dirChangeIdx:N;

        I(backwardPivots) = I(backwardPivots) - B(N - backwardPivots + 1);
    end
end
