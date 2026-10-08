function matrix = createRooftopToTreeBasisTrans(treeBasis, rooftopBasisSize)
    treeBasis = treeBasis(:);
    numRows = numel(treeBasis);

    rowIndices = (1:numRows)'; 
    colIndices = abs(treeBasis); 
    values = sign(treeBasis); 

    matrix = sparse(rowIndices, colIndices, values, numRows, rooftopBasisSize);
end
