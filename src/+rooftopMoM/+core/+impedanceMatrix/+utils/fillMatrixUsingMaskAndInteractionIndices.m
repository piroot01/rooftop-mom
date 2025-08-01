function values = fillMatrixUsingMaskAndInteractionIndices(values, mask, interactionIndices)
% fillMatrixUsingMaskAndInteractionIndices: In-place distribution of interaction values across the matrix using a mask.
%
%   Inputs:
%       values: NxN matrix with initial values at interactionIndices.
%       mask: NxN uint32 mask: non-zero entries index into interaction values.
%       interactionIndices: Linear indices specifying unique values in 'values'.
%
%   Output:
%       values: Fully populated matrix (modified in-place).

    % Extract unique interaction values.
    uv = values(interactionIndices);

    % Identify masked positions.
    nz = (mask ~= 0);

    % Distribute values for non-zero mask entries.
    values(nz)  = uv(mask(nz));

    % Distribute values for zero mask entries.
    values(~nz) = uv(~nz);
end

