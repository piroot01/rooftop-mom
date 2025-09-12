function result = truncate(matrix, countX, countY)
% truncate: Return the top-left submatrix (truncation operator).
%
%   Inputs
%       matrix:M-by-N numeric matrix.
%       countX: Nonnegative integer: number of rows to keep.
%       countY: Nonnegative integer: number of columns to keep.
%
%   Output
%       result: countX-by-countY numeric matrix.

    result = matrix(1, 1:countY);
end
