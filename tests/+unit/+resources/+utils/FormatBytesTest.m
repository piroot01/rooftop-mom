classdef FormatBytesTest < unit.Base
% FormatBytesTest: Tests of rooftopMoM.resources.utils.formatBytes.

    methods (Test)
        function testUnits(testCase)
            testCase.verifyEqual(rooftopMoM.resources.utils.formatBytes(512), '512 B');
            testCase.verifyEqual(rooftopMoM.resources.utils.formatBytes(1536), '1.5 kB');
            testCase.verifyEqual(rooftopMoM.resources.utils.formatBytes(3 * 2^30), '3.0 GB');
            testCase.verifyEqual(rooftopMoM.resources.utils.formatBytes(Inf), 'Inf B');
        end
    end
end
