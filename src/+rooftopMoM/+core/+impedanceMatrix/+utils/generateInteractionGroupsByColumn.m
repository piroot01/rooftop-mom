function groups = generateInteractionGroupsByColumn(interactionIndices, count)
% generateInteractionGroupsByColumn: Groups linear indices into contiguous runs per column of an N×N impedance matrix.
%
%   Inputs:
%       interactions: 1×K vector of linear indices specifying elements to compute.
%       count: Scalar size of the square matrix (N×N).
%
%   Outputs:
%       groups: 1×M cell array; each cell contains a vector of indices that are contiguous (difference of 1) and within the same column.

    % Sort indices and compute corresponding columns.
    idx = sort(interactionIndices(:));
    cols = floor((idx - 1) / count) + 1;

    % Identify boundaries: new column or non-consecutive index.
    boundaries = [true; diff(cols) ~= 0 | diff(idx) ~= 1];

    % Determine group start/end positions.
    starts = find(boundaries);
    ends   = [starts(2:end)-1; numel(idx)];
    sizes  = ends - starts + 1;

    % Split into cell array of contiguous runs.
    groups = mat2cell(idx, sizes);
end
