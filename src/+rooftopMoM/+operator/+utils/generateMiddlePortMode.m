function mode = generateMiddlePortMode(count, useCeil)
% generateMiddlePortMode: Generates a binary vector selecting the middle port.
%
%   Inputs:
%       count: Number of the basis functions.
%       useFloor: (optional) Logical flag specifying which rounding method to use for odd N:
%           true: use ceil(N/2) (default)
%           false: use floor(N/2)
%
%   Output:
%       mode: N-by-1 logical vector with a single true at the computed middle port index.

    % Validate rounding flag, default to use floor if unspecified or empty
    if nargin < 2 || isempty(useCeil)
        useCeil = true;
    end

    % Determine number of ports and initialize mode vector
    mode = false(count, 1);

    % Compute middle index and set that port active
    if useCeil
        idx = ceil(count / 2);
    else
        idx = floor(count / 2);
    end

    mode(idx) = true;
end
