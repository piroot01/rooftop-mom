function mask = generateSingularityMask(indices)
% generateSingularityMask: Build a binary mask of singularity occurrences.
%
%   Input:
%       indices: 1×M cell array, each cell is a 2×K matrix of basis-function indices (rows = polarity 1/2, columns = positions).
%
%   Output:
%       mask: (S×S×4) logical array, where S = sum(K_i). 
%           mask(:,:,n) is the block-diagonal concatenation of per-cell comparisons for polarity pairing n.

    % Define the four polarity-pairings.
    p1 = [1, 1, 2, 2];
    p2 = [1, 2, 1, 2];
    pCount = numel(p1);

    % Compute the number of columns (K) in each cell to size the blocks:
    sizes = cellfun(@(A) size(A,2), indices);
    maskSize = sum(sizes);

    % Preallocate final mask:
    mask = false(maskSize, maskSize, pCount);

    % Precompute block starting offsets along each dimension.
    edges = [0, cumsum(sizes)];

    % For each polarity-pairing, fill in the block-diagonal entries.
    for p = 1:pCount
        % For each cell, compute its little mask and place into mask(:,:,p).
        for n = 1:numel(indices)
            idx = indices{n};

            % Compare the two chosen polarity rows across all positions:
            smallMask = (idx(p1(p), :).' == idx(p2(p), :));

            % Compute the row/col range in the big mask.
            range = (edges(n) + 1):edges(n + 1);

            mask(range, range, p) = smallMask;
        end
    end
end
