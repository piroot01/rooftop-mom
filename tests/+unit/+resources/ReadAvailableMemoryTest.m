classdef ReadAvailableMemoryTest < unit.Base
% ReadAvailableMemoryTest: Tests of rooftopMoM.resources.readAvailableMemory and its parsers.

    methods (Test)
        function testReadsThisMachine(testCase)
            bytes = rooftopMoM.resources.readAvailableMemory();

            testCase.verifyGreaterThan(bytes, 0);
            testCase.verifyTrue(isfinite(bytes));
        end

        function testReadsWindowsAsTheFunctionMemory(testCase)
            % The Java runtime reads the same available memory as the function memory; it may change a little between
            % the two readings.
            testCase.assumeTrue(ispc && usejava('jvm'), 'Only Windows with the Java runtime reads through Java.');
            [~, systemView] = memory();

            testCase.verifyEqual(rooftopMoM.resources.readAvailableMemory(), systemView.PhysicalMemory.Available, ...
                'RelTol', 0.05);
        end

        function testParsesMeminfoOfLinux(testCase)
            text = sprintf(['MemTotal:       32559532 kB\nMemFree:        10826360 kB\n' ...
                'MemAvailable:   19411080 kB\nBuffers:          510956 kB\n']);

            testCase.verifyEqual(rooftopMoM.resources.utils.parseMeminfo(text), 19411080 * 1024);
            testCase.verifyTrue(isnan(rooftopMoM.resources.utils.parseMeminfo('MemTotal: 1 kB')));
        end

        function testParsesVmStatOfMacOS(testCase)
            text = sprintf(['Mach Virtual Memory Statistics: (page size of 16384 bytes)\n' ...
                'Pages free:                               12000.\nPages active:                            500000.\n' ...
                'Pages inactive:                           30000.\nPages speculative:                         1000.\n' ...
                'Pages throttled:                              0.\n']);

            testCase.verifyEqual(rooftopMoM.resources.utils.parseVmStat(text), 16384 * (12000 + 30000 + 1000));
            testCase.verifyTrue(isnan(rooftopMoM.resources.utils.parseVmStat('Pages free: 1.')));
        end
    end
end
