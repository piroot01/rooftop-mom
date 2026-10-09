function bytes = parseVmStat(text)
% parseVmStat: Returns the available memory from the output of the macOS command vm_stat: the free, inactive and
%   speculative pages times the page size given in its first line.
%
%   Input:
%       text: Character vector of the output.
%
%   Output:
%       bytes: Available memory in bytes, NaN if the page size or the free or inactive pages are missing.

    pageSize = readNumber_(text, 'page size of (\d+) bytes');
    free = readNumber_(text, 'Pages free:\s*(\d+)');
    inactive = readNumber_(text, 'Pages inactive:\s*(\d+)');
    speculative = readNumber_(text, 'Pages speculative:\s*(\d+)');

    % Older systems print no speculative pages.
    if isnan(speculative)
        speculative = 0;
    end

    bytes = pageSize * (free + inactive + speculative);
end

function number = readNumber_(text, pattern)
    value = regexp(text, pattern, 'tokens', 'once');
    number = NaN;

    if ~isempty(value)
        number = str2double(value{1});
    end
end
