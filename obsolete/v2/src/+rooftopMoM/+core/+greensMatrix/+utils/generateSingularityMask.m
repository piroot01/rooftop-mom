function mask = generateSingularityMask(indices, rooftopIndex)
% generateSingularityMask: Binary equality mask for polarity pairings.
%
%   Inputs
%       indices: 2-by-N numeric array.
%       rooftopIndex: (optional) scalar or vector of column indices into `indices` to compare against. Default: 1.
%
%   Output
%       mask: logical array of size N-by-M-by-4, where
%           N = number of basis positions (size(indices,2)),
%           M = numel(rooftopIndex),
%           third dimension enumerates the four polarity pairings:
%               pairing 1: indices(1,:) vs indices(1, rooftopIndex)
%               pairing 2: indices(1,:) vs indices(2, rooftopIndex)
%               pairing 3: indices(2,:) vs indices(1, rooftopIndex)
%               pairing 4: indices(2,:) vs indices(2, rooftopIndex)
%           mask(i,j,k) is true when the k-th polarity pairing matches between
%           indices(:, i) and indices(:, rooftopIndex(j)).

    % Define the four polarity pairings (p1 vs p2)
    p1 = [1, 1, 2, 2];
    p2 = [1, 2, 1, 2];
    pCount = numel(p1);

    if nargin < 2 || isempty(rooftopIndex)
        rooftopIndex = 1;
    end

    if rooftopIndex == 1
        % hyperoptimalised

        if indices(2, 1) == 2
            % x-direction
            mask = [1; 2; 0; 1];
        else
            % y-direction
            mask = [1; double(indices(2, 1)); 0; 1];
        end
    else
        mask = zeros(pCount, 1);

        for p = 1:pCount
            %mask(:, :, p) = indices(p1(p), :).' == indices(p2(p), rooftopIndex);
            idx = find(indices(p1(p), :).'  == indices(p2(p), rooftopIndex), 1, 'first');
            if ~isempty(idx)
                mask(p) = idx;
            end
        end
    end
end
