classdef ConsoleSinkTest < unit.Base
    methods (Test)
        function testWritePrintsToCommandWindow(tc)
            sink = logger.sink.ConsoleSink(); %#ok<NASGU>

            captured = evalc('sink.write(sprintf(''Hello Console\\n''))');

            tc.verifyNotEmpty(captured);
            tc.verifyThat(captured, matlab.unittest.constraints.ContainsSubstring('Hello Console'));
        end
    end
end
