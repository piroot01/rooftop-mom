function values = fillMatrixUsingMaskAndInteractionIndices(values, mask, interactionIndices)
% fillMatrixUsingMaskAndInteractionIndices: Vectorized in-place distribution of interaction values across the matrix.
%
%   Inputs:
%       values: N×N matrix, with “source” values sitting at interactionIndices.
%       mask: N×N uint32 mask; nonzero entries index into interaction values.
%       interactionIndices: :inear indices into values that pick out one entry per mask-ID.
%
%   Output:
%       values: N×N matrix with every position filled according to mask.

    % Pull out the mask IDs and the corresponding “source” values.
    maskIDs = mask(interactionIndices);
    sourceValues = values(interactionIndices);

    % Build a lookup so that lookup(k) == the value for mask==k.
    lookup = accumarray(maskIDs, sourceValues);

    % now index into that lookup with the whole mask matrix.
    values = lookup(mask);
end
