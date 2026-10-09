function bytes = parseMeminfo(text)
% parseMeminfo: Returns the available memory from the text of the Linux file /proc/meminfo, the line
%   MemAvailable: 19411080 kB.
%
%   Input:
%       text: Character vector of the file.
%
%   Output:
%       bytes: Available memory in bytes, NaN if the line is missing.

    value = regexp(text, 'MemAvailable:\s*(\d+)\s*kB', 'tokens', 'once');
    bytes = NaN;

    if ~isempty(value)
        bytes = 1024 * str2double(value{1});
    end
end
