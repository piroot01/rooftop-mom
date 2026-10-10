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
%       moments: 7xC moments of the C cell pairs, divided by the squared cell area and ordered as by
%           singularIntegrals.generatePowerMoments.
%       offsets: 2xC distinct nonnegative offsets of the cell pairs, one per cell pair.
%       scales: 1x2 factors of the electric and of the magnetic generators.
%
%   Output:
%       parts: The struct with the moments added.

    if isempty(offsets)
        return;
    end

    [countX, countY] = size(parts.electric);

    % Find the smallest y offset dY, then the number of y offsets up to the largest one, and their column indices in
    % the generators (index = offset + 1).
    first = min(offsets(2, :));
    width = max(offsets(2, :)) - first + 1;
    columns = first + (1:width);

    % The rows give the combinations of M00, M10, M01, M11, N10, N01 and N11. The first one is for the electric
    % generator, the rest for the magnetic ones. The second row is magneticXX at e (both wedges rising or both
    % falling), the third at e - u and the fourth at e + u. Rows five to seven are magneticYY in the same manner; the
    % fifth row uses M00, because N00 = M00, see Definition 2.3.
    weights = [1, 0, 0, 0, 0, 0, 0
        1, -1, -1, 2, 0, 0, 0
        0, 1, 0, -1, 0, 0, 0
        0, 0, 1, -1, 0, 0, 0
        1, 0, 0, 0, -1, -1, 2
        0, 0, 0, 0, 1, 0, -1
        0, 0, 0, 0, 0, 1, -1] .* [scales(1); repmat(scales(2), 6, 1)];

    % Column n of terms holds term n on the countX x width block of the grid: countX offsets along x times width
    % offsets along y, and for each of them 7 terms, each a combination of the moments.
    terms = zeros(countX * width, 7, 'like', 1i);

    % Converts the two offsets into a single linear index of the block: row dX + 1, column dY - first + 1.
    positions = offsets(1, :) + 1 + (offsets(2, :) - first) * countX;
    terms(positions, :) = (weights * moments).';

    % terms(:, n) takes term n out of the 7, and reshape lays it out as the countX x width block of the offsets.
    block = @(n) reshape(terms(:, n), countX, width);

    % The electric generator needs no shift, so it can be added straight away.
    parts.electric(:, columns) = parts.electric(:, columns) + block(1);

    if ~isempty(parts.magneticXX)
        % All interactions rising - falling, M10 - M11, added at e - u.
        rising = block(3);

        % The generator at d_x = 0 also needs the falling term M01 - M11 of the pair at dX = -1, for which there is no
        % data. By (11) it equals the rising term of the pair at dX = 1 (the two wedge pairs are mirror images), so we
        % use this trick and count the row of dX = 1 twice.
        rising(2, :) = 2 * rising(2, :);

        % All interactions falling - rising, M01 - M11, added at e + u.
        falling = block(4);

        % values starts with the second term, M00 - M10 - M01 + 2 M11, added at e.
        values = block(2);

        % Keep the rows dX = 0 to countX - 2, as the generator has countX - 1 rows (the last rooftop is missing). The
        % rising term goes to e - u, so its rows from dX = 1 on are shifted up by one; dX = 0 would land on d_x = -1,
        % which is not needed.
        values = values(1:(end - 1), :) + rising(2:end, :);

        % Similar logic is applied to the falling term: it goes to e + u, so its rows dX = 0 to countX - 3 are shifted
        % down by one, onto d_x = 1 to countX - 2. The last two rows would land outside the generator, and d_x = 0 got
        % its falling term above.
        values(2:end, :) = values(2:end, :) + falling(1:(end - 2), :);

        % And finally add the result to the columns of this chunk.
        parts.magneticXX(:, columns) = parts.magneticXX(:, columns) + values;
    end

    if ~isempty(parts.magneticYY)
        % Here we perform the same as in the magneticXX case, along y with u = e_y, but with two twists: the block holds
        % only the columns of this chunk, and the shifts by one column reach outside of it. For example, countY = 5, so
        % the generator has d_y = 0 to 3, and the last chunk has first = 3 and width = 2 (dY = 3, 4):
        %   column n of values      1    2    3    4
        %   offset first + n - 2    2    3    4    5
        %   block(5) from dY             3    4
        %   rising from dY          3    4
        %   block(7) from dY                  3    4
        %   kept (inside)           yes  yes  no   no
        % Column 1 adds to d_y = 2, which the previous chunk adds to as well.

        % All interactions rising - falling, N10 - N11, added at e - u.
        rising = block(6);

        if first <= 1 && first + width > 1
            % Here we use the same trick (11), but unlike x, the block holds only the columns first to
            % first + width - 1: dY = 1 is its column 2 - first, and only if this chunk contains it.
            rising(:, 2 - first) = 2 * rising(:, 2 - first);
        end

        % Column n of values is the column first + n - 1 of the generator, first - 1 to first + width.
        values = zeros(countX, width + 2, 'like', 1i);

        % The term at e unshifted, the rising term one column left (e - u) and the falling term, N01 - N11, one column
        % right (e + u).
        values(:, 2:(end - 1)) = block(5);
        values(:, 1:(end - 2)) = values(:, 1:(end - 2)) + rising;
        values(:, 3:end) = values(:, 3:end) + block(7);

        % Keep the columns of values whose offset first + n - 2 lies in the generator, 0 to countY - 2: the first chunk
        % drops offset -1, and a chunk at the end of the grid drops the offsets beyond countY - 2.
        inside = max(1, 2 - first):min(width + 2, countY - first);

        % Convert them to the generator columns, index = offset + 1.
        target = first + inside - 1;

        parts.magneticYY(:, target) = parts.magneticYY(:, target) + values(:, inside);
    end
end
