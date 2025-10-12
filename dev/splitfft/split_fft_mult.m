function b = split_fft_mult(pre, x)
% b = toeplitz_matvec(pre, x)
% Input:
%   pre : struct from toeplitz_precompute
%   x   : vector of length prod(N) or array of shape N
% Output:
%   b   : vector same shape (returns column vector)
    N = pre.N;
    d = pre.d;

    xA = reshape(x, N);

    % forward tree: iteratively build even/odd branches for each dimension
    % start with one cell
    cells = {xA};
    for dim = 1:d
        newcells = cell(0);
        Pk = pre.Pgrid{dim};
        for c = 1:numel(cells)
            A = cells{c};
            A_e = fft(A, [], dim);                       % even-branch
            tmp = Pk .* A;                % P .* A
            A_o = fft(tmp, [], dim);                    % odd-branch
            A_o = circshift(A_o, -1, dim);              % align as in 2D code
            newcells{end+1} = A_e; %#ok<AGROW>
            newcells{end+1} = A_o; %#ok<AGROW>
        end
        cells = newcells;
    end

    % NOTE: Up to this point it seems fine.
    %cells{:}

    % multiply each parity block by the corresponding precomputed hc_sub
    parityCount = 2^d;
    for mask = 0:(parityCount-1)
        idx = mask + 1;
        cells{idx} = pre.hc_sub{idx} .* cells{idx};
    end

    % NOTE: Still fine.
    %cells{:}

    % TODO: The part bellow is fucked up.
    % inverse tree: combine pairs back in reverse dimension order
    for dim = d:-1:1
        newcells = cell(0);
        Pk = pre.Pgrid{dim};
        for i = 1:2:numel(cells)
            even = cells{i};
            odd  = cells{i+1};
            he = ifft(even, [], dim);
            ho = ifft(odd,  [], dim);
            combined = he + Pk .* ho;
            newcells{end+1} = combined; %#ok<AGROW>
        end
        cells = newcells;
    end

    % final scaling 1 / 2^d and return real vector
    out = cells{1} / (2^d);
    b = out(:);
end