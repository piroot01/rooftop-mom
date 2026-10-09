function parts = accumulateMoments(parts, moments, offsets, scales)
% accumulateMoments: Adds the moments of cell pairs to the generators of the Green's matrix (the near-field note in
%   docs/nearfield, Generators from moments): scales(1) times M00 to the electric generator and scales(2) times the
%   following combinations to the magnetic ones. A rooftop rises, xi, on its first cell and falls, 1 - xi, on its
%   second, so with u the unit offset along a direction, the cell pair at the offset e adds to the generator of that
%   direction at three offsets:
%       both rising or both falling, at e:                M00 - M10 - M01 + 2 M11;
%       reference rising and the other falling, at e - u: M10 - M11;
%       reference falling and the other rising, at e + u: M01 - M11.
%   At 0 the last term comes from the pair at -u, which the reflection x -> -x maps onto the pair at u with the weight
%   M10 - M11.
%
%   Inputs:
%       parts: Struct of the generators electric (countX x countY), magneticXX and magneticYY, the latter empty for a
%           missing rooftop direction.
%       moments: 7xC moments of the cell pairs, divided by the squared cell area and ordered as by
%           singularIntegrals.generatePowerMoments.
%       offsets: 2xC distinct nonnegative offsets of the cell pairs.
%       scales: 1x2 factors of the electric and of the magnetic generators.
%
%   Output:
%       parts: The struct with the moments added.

    % A generator one cell wide keeps its orientation when indexed, so the values take the shape of the entries.
    index = sub2ind(size(parts.electric), offsets(1, :) + 1, offsets(2, :) + 1);
    entries = parts.electric(index);
    parts.electric(index) = entries + reshape(scales(1) * moments(1, :), size(entries));

    % The y moments take the place of the x moments along y; M00 is common to both.
    names = {'magneticXX', 'magneticYY'};
    rows = {1:4, [1, 5, 6, 7]};

    for d = 1:2
        if isempty(parts.(names{d}))
            continue;
        end

        along = offsets(d, :);
        scaled = scales(2) * moments(rows{d}, :);
        terms = {
            scaled(1, :) - scaled(2, :) - scaled(3, :) + 2 * scaled(4, :), along
            (scaled(2, :) - scaled(4, :)) .* (1 + (along == 1)), along - 1
            scaled(3, :) - scaled(4, :), along + 1
        };

        for n = 1:size(terms, 1)
            [values, target] = terms{n, :};
            inside = target >= 0 & target < size(parts.(names{d}), d);
            position = offsets(:, inside);
            position(d, :) = target(inside);
            index = sub2ind(size(parts.(names{d})), position(1, :) + 1, position(2, :) + 1);
            entries = parts.(names{d})(index);
            parts.(names{d})(index) = entries + reshape(values(inside), size(entries));
        end
    end
end
