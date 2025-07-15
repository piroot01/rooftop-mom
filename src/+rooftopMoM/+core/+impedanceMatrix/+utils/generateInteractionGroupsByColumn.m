function groups = generateInteractionGroupsByColumn(interactions, N)
% generateInteractionGroupsByColumn: Groups interaction indices by column of the impedance matrix.
%
%   Inputs:
%       interactions: 1×K array of linear indices specifying elements to compute.
%       N: Size of the square matrix (N×N).
%
%   Output:
%       groups: 1×M cell array of index vectors, each corresponding to one column.

    sortedInteractions = sort(interactions(:).');

    % Determine column number for each index.
    columnNumbers = floor(double(sortedInteractions - 1) / N) + 1;

    % Find boundaries where the column changes.
    boundaries = [diff(columnNumbers) ~= 0, true];
    boundaryIndices = find(boundaries);

    groupSizes = diff([0, boundaryIndices]);

    % Split sorted indices into groups by column.
    groups = mat2cell(sortedInteractions(:), groupSizes);
end
