function integrals = integrateFacePairs(context, gridSize, shift, weights, staticIntegral)
% integrateFacePairs: Integrates the Green's function exp(-ikR) / R with polynomial weights over the face pairs of a
%   grid. The pair (i, j), counted from 0, is the face (i, j) + shift paired with the face at the origin.
%   On the coincident pair the static part 1 / R is not integrated numerically: the quadrature integrates the regular
%   part (exp(-ikR) - 1) / R and staticIntegral, the analytic integral of the static part, is added.
%
%   Inputs:
%       context: rooftopMoM.assembler.AssemblerContext.
%       gridSize: 1x2 vector [rowCount, columnCount] of the grid.
%       shift: 2x1 face offset added to every grid position.
%       weights: BxB quadrature weights times the weight function, the shifted face along the rows.
%       staticIntegral: Integral of the weight function times 1 / R over the coincident pair.
%
%   Output:
%       integrals: rowCount x columnCount integrals over the unit squares of both faces (without the face areas).

    rooftops = context.rooftops();
    primitive = rooftops.primitive();
    [pointsU, pointsV] = rooftops.quadrature().points();
    wavenumber = context.wavenumber();
    pairCount = numel(weights);
    count = prod(gridSize);

    % The kernel step holds the distances and three complex arrays at once, 56 bytes per point pair.
    chunkLength = context.chunkLength(56 * pairCount);
    chunkCount = ceil(count / chunkLength);

    logger = rooftopMoM.logging.Logger.instance();
    isTrace = logger.isEnabled(rooftopMoM.logging.LogLevel.TRACE);

    integrals = complex(zeros(gridSize));

    for n = 1:chunkCount
        indices = ((n - 1) * chunkLength + 1):min(n * chunkLength, count);
        [i, j] = ind2sub(gridSize, indices);
        offsetX = reshape(i - 1 + shift(1), 1, 1, []);
        offsetY = reshape(j - 1 + shift(2), 1, 1, []);

        % Direct differences keep the distance of coincident points exactly zero.
        distance = sqrt((primitive.sizeX() * (pointsU.' - pointsU + offsetX)).^2 + ...
            (primitive.sizeY() * (pointsV.' - pointsV + offsetY)).^2);
        kernel = exp(-1i * wavenumber * distance) ./ distance;

        coincident = find(offsetX == 0 & offsetY == 0);

        if ~isempty(coincident)
            % expm1 avoids the cancellation in exp(-ikR) - 1 for small kR; -ik is the limit at R = 0.
            regular = expm1(-1i * wavenumber * distance(:, :, coincident)) ./ distance(:, :, coincident);
            regular(distance(:, :, coincident) == 0) = -1i * wavenumber;
            kernel(:, :, coincident) = regular;
        end

        integrals(indices) = weights(:).' * reshape(kernel, pairCount, []);

        if isTrace
            logger.trace('Integrated chunk %d of %d.', n, chunkCount);
        end
    end

    coincidentPosition = 1 - shift(:).';

    if all(coincidentPosition >= 1 & coincidentPosition <= gridSize)
        index = sub2ind(gridSize, coincidentPosition(1), coincidentPosition(2));
        integrals(index) = integrals(index) + staticIntegral;
    end
end
