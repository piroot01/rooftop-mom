classdef ConsoleSink < logger.sink.Sink
    methods
        function write(~, formttedMessage)
            fprintf('%s', formttedMessage);
        end
    end
end
