function index = binarySearchSorted(indices, value)
% binarySearchSorted: Finds the index of the largest element less or equal value.
%
%   Inputs:
%       indices: 1×N sorted array of values.
%       value: Scalar value to search for.
%
%   Output:
%       index: Position of the largest element less or equal value (returns 0 if all elements are greater than value).

    lo = int32(1);
    hi = int32(numel(indices));

    while lo <= hi
        mid = lo + bitshift(hi - lo, -1);

        if indices(mid) <= value
            lo = mid + 1;
        else
            hi = mid - 1;
        end
    end

    index = hi;
end
