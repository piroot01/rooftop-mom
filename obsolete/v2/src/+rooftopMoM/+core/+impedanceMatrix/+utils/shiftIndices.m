function shifted = shiftIndices(indices, offset)
% shiftIndices: Shifts an index mask by a scalar offset.
%
%   Inputs:
%       indices: A numerical array of indices.
%       offset: An integer value to add to each element of the mask.
%
%   Output:
%       shifted: An array of the same size as 'mask' with the shifted indices.

    shifted = indices + offset;
end
