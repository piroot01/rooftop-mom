classdef (Abstract) Sink < handle
% Sink: Destination of the log lines.

    methods (Abstract)
        % write: Writes one line (without the line break).
        write(obj, line);
    end
end
