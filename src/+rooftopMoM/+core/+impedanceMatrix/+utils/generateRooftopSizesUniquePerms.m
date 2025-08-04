function [sizes, indices] = generateRooftopSizesUniquePerms(values, precision)
% generateRooftopSizesUniquePerms: Compute unique rooftop size permutations.
%
%   Inputs:
%       values: 1×N cell array of structs from UniformRooftops.
%       precision: Error of determining the size equality.
%
%   Output:
%       sizes: 2xM array of unique [length; height] pairs (no duplicates).
%       indices: 1xN vector such that sizes(:, indices(i)) corresponds to the size of values{i}.

    % Default to the 10 decimal places.
    if nargin < 2 || isempty(precision)
        precision = 10;
    end

    % Extract sizes into an Nx2 matrix (rows = [length, height])
    allSizes = cell2mat(cellfun(@(rooftop) [rooftop.length; rooftop.height], values, 'UniformOutput', false));

    % Find unique permutations.
    [sizes, ~, indices] = unique(round(allSizes.', precision), 'rows', 'stable');

    % Cast to uint32.
    indices = uint32(indices);
end
