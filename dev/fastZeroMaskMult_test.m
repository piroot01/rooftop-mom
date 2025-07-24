N = 10;
M = 2;

A = rand(N, N);

mask = true(1, N);
missing = randperm(N, M);
mask(missing) = false;

A(~mask, :) = 0;
A(:, ~mask) = 0;

v = rand(N, 1);

B = rand(N, N);

% ref
tic;

Br = B * v;

toc;

% old
tic;

Bo = A * v

matlabsImpl = toc

tic;

ve = zeros(N, 1);
ve(~mask) = v(~mask);

A * ve;

zeroVImpl = toc

% competition
tic;

Bc = A(mask, mask) * v(mask);

maskedImpl = toc

% new version
tic;

Bn = efficient_sparse_mv(A, v, missing);

forImpl = toc

tic;

Bm = masked_matvec_sparse(A, v, mask)

sparseImpl = toc

tic;

mv_bitmask_dot(A, v, mask);

mvBitmaskImpl = toc

function y = masked_matvec_sparse(A, x, mask)
% A    : n×n double
% mask : n×1 logical
% x    : n×1 vector

    n = size(A,1);

    % build a sparse diagonal of 0/1
    D = spdiags(mask, 0, n, n);      % only k nonzeros on the diagonal
    
    % D*A*D is sparse with only k×k nonzero block
    % multiply by x in one go
    y = D * (A * (D * x));
end


function y = efficient_sparse_mv(M, v, active_idx)
% EFFICIENT_SPARSE_MV Computes M*v for a matrix with known zero rows/columns.
%
%   y = EFFICIENT_SPARSE_MV(M, v, active_idx) calculates the matrix-vector
%   product y = M*v efficiently by only using the rows and columns
%   specified in the 'active_idx' vector.
%
%   This function AVOIDS creating intermediate sub-matrices or sub-vectors,
%   thus minimizing memory allocations. It directly calculates the non-zero
%   elements of the result vector 'y'.
%
%   Inputs:
%       M          - The n-by-n source matrix.
%       v          - The n-by-1 source vector.
%       active_idx - A vector containing the 1-based indices of the rows
%                    and columns that are NOT entirely zero.
%
%   Output:
%       y          - The n-by-1 result vector.

% --- Input Validation ---
n = size(M, 1);
if size(M, 2) ~= n || length(v) ~= n
    error('Matrix must be square and dimensions must match vector length.');
end

% --- Algorithm ---

% 1. Pre-allocate the final result vector with zeros.
% This is the only allocation and is necessary for the output.
% The elements corresponding to zero rows will correctly remain zero.
y = zeros(n, 1);

% 2. Loop only over the active rows to calculate the non-zero y values.
% We iterate using an index `k` from 1 to the number of active indices.
for k = 1:length(active_idx)
    % Get the actual row index 'i' from the active set
    i = active_idx(k);

    % This will hold the sum for y(i)
    current_sum = 0;

    % 3. Loop only over the active columns for the dot product.
    % This computes the dot product of row `i` and vector `v`, but skips
    % all terms where the column in M (and element in v) are zero.
    for l = 1:length(active_idx)
        % Get the actual column index 'j' from the active set
        j = active_idx(l);
        
        % Accumulate the product of the active elements
        current_sum = current_sum + M(i, j) * v(j);
    end

    % 4. Assign the calculated sum to the correct position in the result vector.
    y(i) = current_sum;
end

end

function y = mv_bitmask_dot(M, v, ai)
% ai : vector of active (1‑based) indices
    n = size(M,1);
    y = zeros(n,1);
    ai = ai(:);
    subv = v(ai);              % one small allocation
    for k = 1:numel(ai)
        i = ai(k);
        % BLAS‐accelerated dot of length‑k
        y(i) = M(i, ai) * subv;
    end
end

function y = masked_matvec_bitwise(M, v, mask)
% MASKED_MATVEC_BITWISE   Compute y = M*v only on rows/cols where mask==true
%
%   y = MASKED_MATVEC_BITWISE(M, v, mask)
%     M    - n×n matrix
%     v    - n×1 vector
%     mask - n×1 logical (or numeric 0/1) vector indicating active rows/cols
%
%   This uses nested while‑loops to scan through the active indices,
%   avoiding any intermediate submatrix or subvector allocations.

    %-- validate inputs
    [n1,n2] = size(M);
    if n1~=n2 || numel(v)~=n1 || numel(mask)~=n1
        error('M must be n×n and v,mask must have length n');
    end

    mask = logical(mask(:));    % ensure column logical
    y    = zeros(n1,1);         % prezero full output

    rows = mask;                % copy of mask for row scanning
    % outer loop: over each active row i
    while any(rows)
        % pick next active row i
        i = find(rows,1,'first');

        % inner: scan active cols for dot-product
        cols = mask;
        s = 0;
        while any(cols)
            j = find(cols,1,'first');
            s = s + M(i,j) * v(j);
            cols(j) = false;
        end

        y(i)    = s;   % store result
        rows(i) = false;
    end
end

