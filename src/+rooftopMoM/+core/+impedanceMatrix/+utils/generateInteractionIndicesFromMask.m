function indices = generateInteractionIndicesFromMask(mask)
% generateInteractionIndicesFromMask: Finds the first linear index for each unique value in a matrix.
%   This function processes an input matrix 'mask' containing non-negative integer 
%   values (IDs). For each unique ID, it finds the linear index of its first
%   appearance in the matrix when scanned in column-major order.
%
%   Input:
%       mask: An [NxM] matrix of non-negative integer IDs.
%
%   Output:
%       indices: A [Px1] column vector of linear indices.
%           Entries for IDs not present in the mask will be 0.
%

    % Handle empty input gracefully.
    if isempty(mask)
        indices = uint32([]);
        return;
    end

    linearIndices = (1:numel(mask)).';
    indices = accumarray(mask(:), linearIndices, [], @min);

    indices = uint32(indices);
end
