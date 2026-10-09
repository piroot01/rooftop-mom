function text = formatBytes(bytes)
% formatBytes: Returns a number of bytes as a short text in B, kB, MB, GB or TB, powers of 1024.
%
%   Input:
%       bytes: Number of bytes (Inf and NaN allowed).
%
%   Output:
%       text: Character vector, e.g. '1.5 GB'.

    units = {'B', 'kB', 'MB', 'GB', 'TB'};

    if ~isfinite(bytes)
        text = sprintf('%g B', bytes);
        return;
    end

    exponent = min(max(floor(log2(max(bytes, 1)) / 10), 0), numel(units) - 1);

    if exponent == 0
        text = sprintf('%d B', round(bytes));
    else
        text = sprintf('%.1f %s', bytes / 1024^exponent, units{exponent + 1});
    end
end
