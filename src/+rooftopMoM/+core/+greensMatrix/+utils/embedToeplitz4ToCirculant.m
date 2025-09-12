function H = embedToeplitz4ToCirculant(Tgen, M, N)
% embedToeplitz4ToCirculant  Circulant-style embedding into a 2M x 2N matrix.
%
%   H = embedToeplitz4ToCirculant(Tgen, M, N)
%
%   This implements the same embedding as the provided class:
%     - pads/truncates the input into a (M+1) x (N+1) core matrix
%     - uses indices rowsTop = 1:M, rowsBottom = (M+1):-1:2,
%       colsLeft = 1:N, colsRight = (N+1):-1:2
%     - assembles the four quadrants into a 2M x 2N matrix
%
%   The function preserves the data class/complexness of Tgen and
%   truncates if Tgen is larger than (M+1)x(N+1).

    % input sizes
    [M0, N0] = size(Tgen);

    % create core matrix of size (M+1) x (N+1), preserving class
    core = zeros(M+1, N+1, 'like', Tgen);

    % copy (or truncate) into core
    copyM = min(M0, M+1);
    copyN = min(N0, N+1);
    core(1:copyM, 1:copyN) = Tgen(1:copyM, 1:copyN);

    % precompute indices (same as your class)
    rowsTop   = 1:M;
    rowsBottom = (M+1):-1:2;
    colsLeft  = 1:N;
    colsRight = (N+1):-1:2;

    % allocate result preserving class
    H = zeros(2*M, 2*N, 'like', Tgen);

    % fill quadrants
    H(1:M,       1:N      ) = core(rowsTop,   colsLeft);
    H(M+1:end,   1:N      ) = core(rowsBottom,colsLeft);
    H(1:M,       N+1:end  ) = core(rowsTop,   colsRight);
    H(M+1:end,   N+1:end  ) = core(rowsBottom,colsRight);
end
