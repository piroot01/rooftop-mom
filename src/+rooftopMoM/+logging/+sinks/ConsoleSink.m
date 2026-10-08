classdef ConsoleSink < rooftopMoM.logging.sinks.Sink
% ConsoleSink: Writes the log lines to the command window.

    methods
        function write(~, line)
        % write: Writes the line to the command window.
        %
        %   Input:
        %       line: Character vector of the formatted line.

            fprintf(1, '%s\n', line);
        end
    end
end
