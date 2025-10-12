function pre = split_fft_mult_prep(ht)
% pre = toeplitz_precompute(ht)
% Input:
%   ht : d-dimensional Toeplitz kernel array with size N = [N1 N2 ... Nd]
% Output:
%   pre : struct with fields
%       - N      : original kernel sizes
%       - d      : dimensionality
%       - hc_sub : cell array of length 2^d with subsampled spectral blocks
%                 each block has size N (i.e. size(ht))
%       - Pgrid  : cell array of d arrays sized like ht containing the phase factors
%       - circSize : vector of circulant sizes (2*Ni)
%
    N = size(ht);
    d = numel(N);

    % circulant sizes (use 2*Ni as in your second method)
    circSize = 2 * N;

    % Build the d-dimensional extended circulant array hc (zeros except mirrored blocks)
    % To do that we create index grids and map each coordinate to kernel index or "hole" (Ni+1)
    idxCells = cell(1,d);
    for k = 1:d
        idxCells{k} = 1:circSize(k);        % positions in the circulant dimension
    end
    [G{1:d}] = ndgrid(idxCells{:});


    % map each G{k} to an index in ht or to 0 if it's the center (Ni+1)
    mapped = cell(1,d);
    validMask = true(size(G{1}));
    for k = 1:d
        Nk = N(k);
        J = G{k};
        Mapped = zeros(size(J));
        % j <= Nk -> mapped = j
        % j == Nk+1 -> hole -> leave 0
        % j >= Nk+2 -> mapped = 2*Nk - j + 2
        leMask = (J <= Nk);
        eqMask = (J == Nk + 1);
        gtMask = (J >= Nk + 2);
        Mapped(leMask) = J(leMask);
        Mapped(eqMask) = 0;                    % hole
        Mapped(gtMask) = 2*Nk - J(gtMask) + 2; % mirrored indices
        mapped{k} = Mapped;
        validMask = validMask & (Mapped > 0);
    end

    % Prepare hc and fill only valid positions from ht
    hc = zeros(circSize);
    % find linear indices of valid positions
    validIdx = find(validMask);
    if ~isempty(validIdx)
        % get the subscript indices in the circulant grid for validIdx
        % TODO: Potentially can be faster.
        %[subs{1:d}] = ind2sub(circSize, validIdx);
        % build corresponding indexes into ht
        ht_subs = cell(1,d);
        for k = 1:d
            % mapped{k} is an array same size as hc
            ht_subs{k} = mapped{k}(validIdx);
        end
        % build linear indices into ht and assign
        % TODO: May be done more efficiently.
        ht_lin = sub2ind(N, ht_subs{:});
        % NOTE: The embedding seems fine.
        hc(validIdx) = ht(ht_lin);
    end

    % spectral transform of hc
    hc_fft = fftn(hc);

    % extract the 2^d parity blocks: for each parity vector p in {0,1}^d
    % indices along each dim: start = 1 + p (i.e. 1 or 2), step = 2 => start:2:end
    parityCount = 2^d;
    hc_sub = cell(parityCount,1);
    for mask = 0:(parityCount-1)
        % NOTE: Gives the dim in a binary manner.
        bits = bitget(mask, d:(-1):1); % bitget gives LSB first, corresponds to dims 1..d
        slices = cell(1,d);
        for k = 1:d
            start = 1 + bits(k);
            slices{k} = start:2:circSize(k);
        end
        % use dynamic indexing to extract
        hc_sub{mask+1} = hc_fft(slices{:});
    end

    % precompute phase factors P_k on original N(k) grid, then expand to full grid shape
    Pgrid = cell(1,d);
    for k = 1:d
        Nk = N(k);
        vec = exp(1i * pi * (0:(Nk-1)).' / Nk);    % column vector
        % create shape to reshape into proper dims
        shape = ones(1,d);
        shape(k) = Nk;
        reshaped = reshape(vec, shape);
        % replicate to N
        reps = N ./ shape;
        Pgrid{k} = repmat(reshaped, reps);
    end

    pre.N = N;
    pre.d = d;
    pre.hc_sub = hc_sub;
    pre.Pgrid = Pgrid;
    pre.circSize = circSize;
end