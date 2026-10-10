function parts = integrateLattice_(parts, context, scales, rule, counts, excluded)
% integrateLattice_: Integrates the moments of the far cell pairs, the offsets 0 to counts - 1 except those below
%   excluded in both directions, with a rule of one part per direction, generateReducedRule(order, [1, 1]), and adds
%   them to the generators.
%   The rule has the nodes nu = p - 1 on the left half and nu = p on the right half, with p the order Gauss-Legendre
%   points on [0, 1]. At the offset dX the kernel needs x = a (nu - dX), which is a (p - m) with m = dX + 1 on the left
%   half and m = dX on the right half. The right half of the offset dX + 1 thus uses the same points as the left half
%   of the offset dX:
%       offset dX       0        1        2        3
%       left half       m = 1    m = 2    m = 3    m = 4
%       right half      m = 0    m = 1    m = 2    m = 3
%   Along y likewise. All offsets together need the kernel only on the lattice of the blocks m = 0, ..., countX times
%   n = 0, ..., countY, order x order points each: order^2 (countX + 1)(countY + 1) kernel values instead of the
%   (2 order)^2 countX countY of integrating every offset on its own, about 4x fewer. The values themselves are the
%   same, so are the moments. The lattice is processed in tiles of offsets, one chunk of memory each.

    primitive = context.rooftops().primitive();
    sizes = [primitive.sizeX(), primitive.sizeY()];
    wavenumber = context.wavenumber();
    order = numel(rule.nodesX) / 2;

    % The nodes of the right half, nu = p; the left half is the same minus 1.
    points = rule.nodesX(order + (1:order));
    quadrants = splitQuadrants_(rule, order);

    [rowCount, columnCount] = chooseTile_(context.chunkSize(), estimateBlockBytes_(order), counts);
    tileCounts = ceil(counts ./ [rowCount, columnCount]);
    tileCount = prod(tileCounts);

    logger = rooftopMoM.logging.Logger.instance();
    isTrace = logger.isEnabled(rooftopMoM.logging.LogLevel.TRACE);

    for tile = 1:tileCount
        % The offsets of the tile: dX in rangeX and dY in rangeY, tiles running down the grid first.
        [tileX, tileY] = ind2sub(tileCounts, tile);
        rangeX = ((tileX - 1) * rowCount):(min(tileX * rowCount, counts(1)) - 1);
        rangeY = ((tileY - 1) * columnCount):(min(tileY * columnCount, counts(2)) - 1);
        moments = integrateTile_(quadrants, points, sizes, wavenumber, rangeX, rangeY);

        % The near offsets lie on the lattice too, since they share its blocks, but their moments come from the near
        % pass, so they are dropped here. p lies inside (0, 1), so no lattice point has R = 0, not even for the
        % touching offsets.
        [dX, dY] = ndgrid(rangeX, rangeY);
        offsets = [dX(:).'; dY(:).'];
        kept = offsets(1, :) >= excluded(1) | offsets(2, :) >= excluded(2);
        parts = rooftopMoM.assembler.greensMatrix.utils.accumulateMoments(parts, moments(:, kept), ...
            offsets(:, kept), scales);

        if isTrace
            logger.trace('Integrated chunk %d of %d.', tile, tileCount);
        end
    end
end

function moments = integrateTile_(quadrants, points, sizes, wavenumber, rangeX, rangeY)
% integrateTile_: Returns the 7 x (numel(rangeX) * numel(rangeY)) moments of the offsets rangeX x rangeY, dX running
%   fastest, from the kernel on the blocks rangeX(1), ..., rangeX(end) + 1 times rangeY(1), ..., rangeY(end) + 1.

    order = numel(points);
    blocksX = [rangeX, rangeX(end) + 1];
    blocksY = [rangeY, rangeY(end) + 1];

    % The kernel is wanted as an order x order x numel(blocksX) x numel(blocksY) array: the point in x, the point in y,
    % the block in x and the block in y. points - blocksX is the order x numel(blocksX) matrix of the x coordinates
    % a (p - m); reshape moves its columns to the dimension 3, and the y coordinates go to the dimensions 2 and 4. hypot
    % then expands both into the whole array, without forming the grids of coordinates first.
    x = reshape(sizes(1) * (points - blocksX), order, 1, numel(blocksX));
    y = reshape(sizes(2) * (points - blocksY), 1, order, 1, numel(blocksY));
    distance = hypot(x, y);
    kernel = exp(-1i * wavenumber * distance) ./ distance;
    clear distance

    % Every column is now one block, its order^2 points with x running fastest, like the columns of the quadrant
    % weights. reshape only relabels the array, it copies nothing.
    kernel = reshape(kernel, order^2, []);

    % The moments of an offset are the sums of its four quadrants, each the quadrant weights times one block. A left
    % half takes the block of the offset plus 1 (shift 1), a right half the block of the offset itself (shift 0); the
    % offset at the tile position k has its own block at the position k.
    rows = 1:numel(rangeX);
    columns = 1:numel(rangeY);
    moments = zeros(7, numel(rangeX), numel(rangeY), 'like', 1i);

    for quadrant = quadrants
        sums = reshape(quadrant.weights * kernel, 7, numel(blocksX), numel(blocksY));
        moments = moments + sums(:, rows + quadrant.shift(1), columns + quadrant.shift(2));
    end

    moments = reshape(moments, 7, []);
end

function quadrants = splitQuadrants_(rule, order)
% splitQuadrants_: Returns the four quadrants of the rule, its left or right half in x times its left or right half in
%   y, as a 1x4 struct array with the fields:
%       weights: 7 x order^2 weights of the quadrant's nodes, x running fastest like the points of a block.
%       shift: 1x2 shift of the quadrant's block from the offset in x and in y, 1 for a left half, 0 for a right one.

    % The columns of rule.weights run over the nodes (nu_k, eta_l) with k fastest, so this reshape gives
    % weights(:, k, l), the weight of the node (nu_k, eta_l). The first order nodes of each axis are the left half.
    weights = reshape(rule.weights, 7, 2 * order, 2 * order);
    halves = {1:order, order + (1:order)};
    quadrants = struct('weights', cell(1, 4), 'shift', cell(1, 4));

    for halfY = 1:2
        for halfX = 1:2
            quadrant = halfX + 2 * (halfY - 1);
            quadrants(quadrant).weights = reshape(weights(:, halves{halfX}, halves{halfY}), 7, []);
            quadrants(quadrant).shift = [halfX == 1, halfY == 1];
        end
    end
end

function [rowCount, columnCount] = chooseTile_(chunkSize, blockBytes, counts)
% chooseTile_: Returns the numbers of offsets along x and along y of a tile. A tile of r x c offsets holds the
%   (r + 1) x (c + 1) blocks they use. It takes whole columns of the grid while the blocks of two of them fit into the
%   chunk, otherwise part of one column, at least one offset.

    columnBytes = (counts(1) + 1) * blockBytes;

    if 2 * columnBytes <= chunkSize
        rowCount = counts(1);
        columnCount = min(counts(2), floor(chunkSize / columnBytes) - 1);
    else
        rowCount = min(counts(1), max(1, floor(chunkSize / (2 * blockBytes)) - 1));
        columnCount = 1;
    end
end
