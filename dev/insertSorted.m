function indices = insertSorted(originalIndices, n)
% insertSorted: Inserts a value into a sorted array, preserving order.
%
%   Inputs:
%       originalIndices: 1×N sorted array of values.
%       n: Scalar value to insert.
%
%   Output:
%       indices: 1×(N+1) sorted array with n inserted.

    k = binarySearchSorted(originalIndices, n);

    % Insert n into the sorted array
    if k == 0
        indices = [n, originalIndices];
    elseif k == numel(originalIndices)
        indices = [originalIndices, n];
    else
        indices = [originalIndices(1:k), n, originalIndices((k + 1):end)];
    end
end
