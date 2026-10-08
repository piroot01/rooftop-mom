function I = applyKInvMatrix(Q, treePatches)
% applyKInvMatrix: Solves a system of linear equations written in the matrix form K * I = Q.
%   The equation can be solved by applying the inverse of K to the right-hand side: I = K^(-1) * Q.
%   By leveraging the known structure of the matrix, we can accelerate the application of K^(-1)
%   substantially, achieving linear time complexity and essentially constant memory.
%
%   Inputs:
%       Q: Nx1 vector defining the right-hand side of the equation.
%       treePatches: Nx1 or 1xN vector defining how the tree spans the patches, and in what order.
%
%   Output:
%       I: Nx1 vector that solves the equation for the given Q and the K matrix defined by the tree.
%
%   Notes:
%       The K matrix is defined as a Gram matrix where each element is an inner product between a
%       pulse function and the divergence of a rooftop basis function within the tree spanning the
%       plate. Hence, K_{mn} = < P_{m} | div(R^{(T)}_{n}) >. Since the divergence of the rooftop
%       is a Haar wavelet function (a superposition of two pulse functions with opposite polarities),
%       we can predict where a [-1, +1] pair will be placed based on the tree patches.
%       An exception occurs when the pivot (the position of +1 within the row) is larger
%       than the actual matrix size. In this case, only the -1 from the pair is present, requiring us
%       to shift our strategy from forward substitution to backward substitution.
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

    % The first element of I will be the same as Q because the first rooftop within the tree
    % does not have any predecessor.
    I(1) = Q(1);

    if N == 1
        return;
    end

    % Pick the relevant pivot indices from the treePatches.
    pivots = treePatches(2:(end - 1));

    % Find the first occurrence where a pivot index exceeds N.
    idx = find(pivots > N, 1);

    % Determine the breaking points for forward and backward substitution upfront.
    if isempty(idx)
        % No missing row; forward substitution covers the whole array.
        idx = N;
    end

    % Index marking the first pivot of the backwards substitution element.
    dirChangeIdx = idx + 1;

    % Forward Substitution:
    % Compute the usual forward substitution up to idx.
    if idx > 1
        I(2:idx) = Q(1) + cumsum(Q(pivots(1:(idx - 1))));
    end

    % Backward Substitution:
    % Compute only if a direction change occurred within the matrix bounds.
    if dirChangeIdx <= N
        % We proceed with the backward substitution by fixing the last value of I.
        % The multiplication by -1 is necessary because only the -1 from the [-1, 1] pair 
        % is present in the matrix.
        I(N) = (-1) * Q(dirChangeIdx);

        % Compute only if the backward substitution has nonzero steps.
        if N > dirChangeIdx
            I((N - 1):(-1):dirChangeIdx) = I(N) - cumsum(Q((dirChangeIdx + 1):N));
        end
    end
end
