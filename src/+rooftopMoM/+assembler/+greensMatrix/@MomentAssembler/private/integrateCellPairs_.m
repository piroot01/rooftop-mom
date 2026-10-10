function parts = integrateCellPairs_(parts, context, scales, rule, extent, excluded)
% integrateCellPairs_: Integrates the moments of the cell pairs at the offsets 0 to extent - 1, except those below
%   excluded in both directions, with the rule and adds them to the generators, one chunk of offsets at a time.

    primitive = context.rooftops().primitive();
    sizes = [primitive.sizeX(), primitive.sizeY()];
    count = prod(extent);

    % A chunk peaks with the distances and the kernel values of integrateMoments, both Px x Py per offset: the distance
    % is a real double (8 bytes) and the kernel a complex one (16 bytes), 24 bytes per node together. chunkLength is
    % thus the number of offsets whose two arrays fit into the memory budget of the context (chunkSize of Problem).
    % The complex temporaries inside the kernel expression are not counted, so the true peak is briefly higher.
    chunkLength = context.chunkLength(24 * numel(rule.nodesX) * numel(rule.nodesY));
    chunkCount = ceil(count / chunkLength);

    logger = rooftopMoM.logging.Logger.instance();
    isTrace = logger.isEnabled(rooftopMoM.logging.LogLevel.TRACE);

    for n = 1:chunkCount
        % The offsets of the other cell in chunk n, without the excluded ones (the near offsets in the far pass).
        range = ((n - 1) * chunkLength + 1):min(n * chunkLength, count);
        offsets = selectOffsets_(extent, range, excluded);

        % The moments of a chunk are a temporary, so they are freed before the next chunk is integrated.
        if ~isempty(offsets)
            parts = rooftopMoM.assembler.greensMatrix.utils.accumulateMoments(parts, ...
                rooftopMoM.assembler.greensMatrix.utils.integrateMoments(rule, offsets, sizes, context.wavenumber()), ...
                offsets, scales);
        end

        if isTrace
            logger.trace('Integrated chunk %d of %d.', n, chunkCount);
        end
    end
end

function offsets = selectOffsets_(extent, range, excluded)
% selectOffsets_: Returns the 2xC offsets with the linear indices range among the offsets 0 to extent - 1, except
%   those below excluded in both directions, one per cell pair. C is numel(range), or fewer if some are excluded.

    [m, n] = ind2sub(extent, range);
    offsets = [m - 1; n - 1];
    offsets = offsets(:, offsets(1, :) >= excluded(1) | offsets(2, :) >= excluded(2));
end
