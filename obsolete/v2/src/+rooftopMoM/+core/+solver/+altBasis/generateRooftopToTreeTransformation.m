function matrix = generateRooftopToTreeTransformation(treeConverted, rooftopBasisSize)
% generateRooftopToTreeTransformation: Generate a sparse matrix describing the transformation
%   from the rooftop basis into the tree basis.
%
%   Inputs:
%       treeConverted: Representation of the tree in the means of the rooftop basis.
%       rooftopBasisSize: Size of the rooftop basis.
%
%   Output:
%       matrix: Sparse matrix representation of the transformation matrix, filled with O(N)
%           elements.

    treeBasis = treeConverted(:);
    numRows = numel(treeBasis);

    rowIndices = (1:numRows)'; 
    colIndices = abs(treeBasis); 
    values = sign(treeBasis); 

    matrix = sparse(rowIndices, colIndices, values, numRows, rooftopBasisSize);
end

