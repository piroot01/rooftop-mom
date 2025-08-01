function shifted = shiftMask(mask, offset)
% shiftMask: Shifts an index mask by a scalar offset.
%
%   Inputs:
%       mask: A numerical array of indices.
%       offset: An integer value to add to each element of the mask.
%
%   Output:
%       shifted: An array of the same size as 'mask' with the shifted indices.

    shifted = mask + offset;
end
