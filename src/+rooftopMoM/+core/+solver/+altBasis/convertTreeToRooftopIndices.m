function treeConverted = convertTreeToRooftopIndices(tree, countX, countY)
% convertTreeToRooftopIndices: Return the tree in the means of the rooftop
%   basis indices while covering the rooftops orientation.
%
%   Inputs:
%       tree: 1xN matrix describing the tree.
%       countX: Positive integer specifying the number of patches in the x direction.
%       countY: Positive integer specifying the number of patches in the y direction.
%
%   Output:
%       treeConverted: 1xN matrix describing the tree basis functions in terms of
%           rooftop indices.

    dBasis = diff(tree); 

    % 'sourceIndex' corresponds to all elements except the last
    sourceIndices = tree(1:end-1);

    treeConverted = zeros(size(dBasis));

    maskPlus1  = (dBasis == 1);
    maskMinus1 = (dBasis == -1);
    maskOther  = ~(maskPlus1 | maskMinus1); % Everything else

    idxP = sourceIndices(maskPlus1);
    treeConverted(maskPlus1) = idxP - floor(idxP / countX);

    idxM = sourceIndices(maskMinus1);
    treeConverted(maskMinus1) = -(idxM - 1 - floor((idxM - 1) / countX));

    treeConverted(maskOther) = sourceIndices(maskOther) + (countX - 1) * countY;
end
