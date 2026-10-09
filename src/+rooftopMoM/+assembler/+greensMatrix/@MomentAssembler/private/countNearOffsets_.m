function nearCounts = countNearOffsets_(counts, sizes, gap)
% countNearOffsets_: Returns the 1x2 numbers of the near cell offsets along x and along y, those whose cells are less
%   than gap longer cell sides apart. Along each direction they are the offsets 0 to nearCounts - 1.

    nearCounts = zeros(1, 2);

    for d = 1:2
        gaps = max((0:(counts(d) - 1)) - 1, 0) * sizes(d);
        nearCounts(d) = nnz(gaps < gap * max(sizes));
    end
end
