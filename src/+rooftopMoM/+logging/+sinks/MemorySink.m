classdef MemorySink < rooftopMoM.logging.sinks.Sink
% MemorySink: Keeps the log lines in memory, e.g. for tests or to inspect a run afterwards.

    properties (Access = private)
        lines_ cell = {};
    end

    methods
        function write(obj, line)
        % write: Appends the line to the stored lines.
        %
        %   Input:
        %       line: Character vector of the formatted line.

            obj.lines_{end + 1, 1} = line;
        end

        function lines = lines(obj)
            lines = obj.lines_;
        end
    end
end
