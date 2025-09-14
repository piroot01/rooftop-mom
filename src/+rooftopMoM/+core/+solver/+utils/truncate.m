function result = truncate(matrix, sizeX, sizeY)
% truncate: Return the top-left submatrix (truncation operator).
%
%   Inputs
%       matrix:M-by-N numeric matrix.
%       sizeX: Nonnegative integer: number of rows to keep.
%       sizeY: Nonnegative integer: number of columns to keep.
%
%   Output
%       result: countX-by-countY numeric matrix.

    result = matrix(1:sizeX, 1:sizeY);
end
