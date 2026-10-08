classdef TestCaptureSink < logger.sink.Sink
% TestCaptureSink: Simple in-memory sink for unit tests.

    properties
        messages cell = {};
    end

    methods
        function write(obj, formattedMessage)
            % Append the message (keeps original formatting).
            obj.messages{end + 1} = formattedMessage;
        end
    end
end
