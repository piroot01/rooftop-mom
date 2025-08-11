function H = embedToeplitz4ToCirculant(Tgen, M, N)
% embedToeplitz4ToCirculant: Circulant embedding of a 4th-order Toeplitz generator.
%
%   Inputs:
%       Tgen: Generating MxN block (or smaller; it will be zero-padded). 
%           Can be real or complex. Size may be (M-1)xN, Mx(N-1), etc.
%       M, N: Target sizes for the top-left generator block.
%
%   Output:
%       H: (2M - 1) x (2N - 1) circulant-embedded matrix constructed by reflecting Tgen:
%           H(i,j) = Tgen(src_row(i), src_col(j)), where
%           src_row = [1:M, M-1:-1:1], 
%           src_col = [1:N, N-1:-1:1].
%
%   Notes:
%       If Tgen is smaller than MxN it is zero-padded to MxN.
%       If Tgen is larger than MxN it is truncated.

    % Sizes of provided Tgen
    [M0, N0] = size(Tgen);

    % Prepare padded/truncated generator block A of size M x N.
    % Use 'like' to preserve class (and complex-ness) of Tgen without creating double.
    if M0 == M && N0 == N
        % No padding/truncation needed.
        A = Tgen;
    else
        % Allocate minimal pad array.
        A = zeros(M, N, 'like', Tgen);

        Mmin = min(M0, M); 
        Nmin = min(N0, N);

        % Copy available entries.
        A(1:Mmin, 1:Nmin) = Tgen(1:Mmin, 1:Nmin);
    end

    % Build source index vectors that reflect the rows and cols:
    % For rows: [1,2,...,M, M-1, ..., 1] (length 2M - 1)
    src_row = [1:M, (M - 1):(-1):1];

    % For cols: [1,2,...,N, N-1, ..., 1] (length 2N-1)
    src_col = [1:N, (N - 1):(-1):1];

    % Single indexed extraction to form H.
    H = A(src_row, src_col);
end
