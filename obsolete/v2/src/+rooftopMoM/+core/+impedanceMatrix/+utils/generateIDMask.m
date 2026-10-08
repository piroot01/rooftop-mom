function mask = generateIDMask(IDMask)
% generateIDMask: Generate an ID‑match mask weighted by the matching ID.
%
%   Input:
%       IDMask: N×1 vector of uint32 basis function IDs.
%
%   Output:
%       mask: N×N uint32 matrix where
%           mask(i,j) = IDMask(j) if IDMask(i)==IDMask(j)
%                       0         otherwise

    IDMask = IDMask(:);
    matches = (IDMask == IDMask.');
    mask = uint32(matches) .* uint32(IDMask.');
end
