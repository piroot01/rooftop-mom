classdef FaultySink < logger.sink.Sink
% FaultySink: Intentionally throws on write to test logger robustness.

    methods
        function write(~, ~)
            error('unit.helpers.logger.sink.FaultySink: Intentional sink failure for tests.');
        end
    end
end
