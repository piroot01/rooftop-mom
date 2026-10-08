function rooftopIndices = mapTreeBasisToRooftop(treeBasis, countX, countY)
    % Map tree basis to rooftop indices and their orientations

    dBasis = diff(treeBasis); 

    % 'sourceIndex' corresponds to all elements except the last
    sourceIndices = treeBasis(1:end-1);

    rooftopIndices = zeros(size(dBasis));

    maskPlus1  = (dBasis == 1);
    maskMinus1 = (dBasis == -1);
    maskOther  = ~(maskPlus1 | maskMinus1); % Everything else

    idxP = sourceIndices(maskPlus1);
    rooftopIndices(maskPlus1) = idxP - floor(idxP / countX);

    idxM = sourceIndices(maskMinus1);
    rooftopIndices(maskMinus1) = -(idxM - 1 - floor((idxM - 1) / countX));

    rooftopIndices(maskOther) = sourceIndices(maskOther) + (countX - 1) * countY;
end
