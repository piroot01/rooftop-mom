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
%   The terms are laid out on the grid columns between the smallest and the largest offset along y, zero where no
%   offset is given, and added to the generators as blocks of whole columns, shifted by one row for e - u and e + u
%   along x and by one column along y.
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

    if isempty(offsets)
        return;
    end

    [countX, countY] = size(parts.electric);
    first = min(offsets(2, :));
    width = max(offsets(2, :)) - first + 1;
    columns = first + (1:width);

    % The rows give the electric term, then e, e - u and e + u along x and along y, from M00, M10, M01, M11, N10, N01
    % and N11; the y moments take the place of the x moments along y.
    weights = [1, 0, 0, 0, 0, 0, 0
        1, -1, -1, 2, 0, 0, 0
        0, 1, 0, -1, 0, 0, 0
        0, 0, 1, -1, 0, 0, 0
        1, 0, 0, 0, -1, -1, 2
        0, 0, 0, 0, 1, 0, -1
        0, 0, 0, 0, 0, 1, -1] .* [scales(1); repmat(scales(2), 6, 1)];

    % Column n of terms holds term n on the countX x width block of the grid.
    terms = zeros(countX * width, 7, 'like', 1i);
    terms(offsets(1, :) + 1 + (offsets(2, :) - first) * countX, :) = (weights * moments).';
    block = @(n) reshape(terms(:, n), countX, width);

    parts.electric(:, columns) = parts.electric(:, columns) + block(1);

    if ~isempty(parts.magneticXX)
        rising = block(3);
        rising(2, :) = 2 * rising(2, :);
        falling = block(4);
        values = block(2);
        values = values(1:(end - 1), :) + rising(2:end, :);
        values(2:end, :) = values(2:end, :) + falling(1:(end - 2), :);
        parts.magneticXX(:, columns) = parts.magneticXX(:, columns) + values;
    end

    if ~isempty(parts.magneticYY)
        rising = block(6);

        if first <= 1 && first + width > 1
            rising(:, 2 - first) = 2 * rising(:, 2 - first);
        end

        % Column n of values is the column first + n - 1 of the generator, first - 1 to first + width.
        values = zeros(countX, width + 2, 'like', 1i);
        values(:, 2:(end - 1)) = block(5);
        values(:, 1:(end - 2)) = values(:, 1:(end - 2)) + rising;
        values(:, 3:end) = values(:, 3:end) + block(7);
        inside = max(1, 2 - first):min(width + 2, countY - first);
        target = first + inside - 1;
        parts.magneticYY(:, target) = parts.magneticYY(:, target) + values(:, inside);
    end
end
