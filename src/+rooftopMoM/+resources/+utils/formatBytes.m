function text = formatBytes(bytes)
% formatBytes: Returns a number of bytes as a short text in B, kB, MB, GB or TB, powers of 1024. Inf, a missing
%   ceiling, reads 'unlimited' and NaN, an unknown free memory, reads 'unknown'.
%
%   Input:
%       bytes: Number of bytes, Inf or NaN.
%
%   Output:
%       text: Character vector, e.g. '1.5 GB'.

    units = {'B', 'kB', 'MB', 'GB', 'TB'};

    if isinf(bytes)
        text = 'unlimited';
        return;
    elseif isnan(bytes)
        text = 'unknown';
        return;
    end

    exponent = min(max(floor(log2(max(bytes, 1)) / 10), 0), numel(units) - 1);

    if exponent == 0
        text = sprintf('%d B', round(bytes));
    else
        text = sprintf('%.1f %s', bytes / 1024^exponent, units{exponent + 1});
    end
end
