classdef LoggerTest < unit.Base
% LoggerTest: Tests of rooftopMoM.logging.Logger.

    methods (Test)
        function testInstanceIsShared(testCase)
            testCase.verifySameHandle(rooftopMoM.logging.Logger.instance(), rooftopMoM.logging.Logger.instance());
        end

        function testDefaultLevelIsWarning(testCase)
            testCase.verifyEqual(rooftopMoM.logging.Logger.instance().level(), rooftopMoM.logging.LogLevel.WARNING);
        end

        function testWritesMessagesUpToTheLevel(testCase)
            sink = testCase.captureLines_(rooftopMoM.logging.LogLevel.INFO);
            logger = rooftopMoM.logging.Logger.instance();

            logger.error('Error.');
            logger.warning('Warning.');
            logger.info('Info.');
            logger.debug('Debug.');
            logger.trace('Trace.');

            lines = sink.lines();
            testCase.verifyNumElements(lines, 3);
            testCase.verifySubstring(lines{3}, '[INFO]');
        end

        function testOffWritesNothing(testCase)
            sink = testCase.captureLines_(rooftopMoM.logging.LogLevel.OFF);
            rooftopMoM.logging.Logger.instance().error('Error.');
            testCase.verifyEmpty(sink.lines());
        end

        function testIsEnabled(testCase)
            testCase.captureLines_(rooftopMoM.logging.LogLevel.DEBUG);
            logger = rooftopMoM.logging.Logger.instance();

            testCase.verifyTrue(logger.isEnabled(rooftopMoM.logging.LogLevel.DEBUG));
            testCase.verifyFalse(logger.isEnabled(rooftopMoM.logging.LogLevel.TRACE));
        end

        function testLineFormat(testCase)
            sink = testCase.captureLines_(rooftopMoM.logging.LogLevel.DEBUG);
            rooftopMoM.logging.Logger.instance().debug('Converged in %d iterations.', 12);

            lines = sink.lines();
            testCase.verifyMatches(lines{1}, ['^\[\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\.\d{3}\] \[DEBUG\] ' ...
                '\[LoggerTest\.testLineFormat\]: Converged in 12 iterations\.$']);
        end

        function testTimestampIsTheCurrentTime(testCase)
            sink = testCase.captureLines_(rooftopMoM.logging.LogLevel.INFO);
            rooftopMoM.logging.Logger.instance().info('Now.');

            lines = sink.lines();
            written = datetime(lines{1}(2:24), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');
            testCase.verifyLessThan(abs(seconds(datetime('now') - written)), 1);
        end

        function testCallerOfLocalFunction(testCase)
            sink = testCase.captureLines_(rooftopMoM.logging.LogLevel.WARNING);
            warnFromLocalFunction_();

            lines = sink.lines();
            testCase.verifySubstring(lines{1}, '[LoggerTest>warnFromLocalFunction_]: ');
        end

        function testCallerOfMethodInClassFolder(testCase)
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture()).Folder;
            classFolder = fullfile(folder, '@CallerProbe');
            mkdir(classFolder);
            writelines({'classdef CallerProbe', '    methods', '        warn(obj);', '    end', 'end'}, ...
                fullfile(classFolder, 'CallerProbe.m'));
            writelines({'function warn(~)', '    rooftopMoM.logging.Logger.instance().warning(''Hello.'');', 'end'}, ...
                fullfile(classFolder, 'warn.m'));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(folder));
            sink = testCase.captureLines_(rooftopMoM.logging.LogLevel.WARNING);

            probe = CallerProbe();
            probe.warn();

            lines = sink.lines();
            testCase.verifySubstring(lines{1}, '[CallerProbe.warn]: Hello.');
        end

        function testEverySinkReceivesTheLine(testCase)
            first = rooftopMoM.logging.sinks.MemorySink();
            second = rooftopMoM.logging.sinks.MemorySink();
            logger = rooftopMoM.logging.Logger.instance();
            logger.configure(rooftopMoM.logging.LogLevel.INFO, {first, second});

            logger.info('Both.');

            testCase.verifyEqual(first.lines(), second.lines());
            testCase.verifyNumElements(first.lines(), 1);
        end

        function testFailingSinkRaisesWarning(testCase)
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture()).Folder;
            logger = rooftopMoM.logging.Logger.instance();

            % A folder cannot be opened as a file, so every write fails.
            logger.configure(rooftopMoM.logging.LogLevel.WARNING, {rooftopMoM.logging.sinks.FileSink(folder)});

            testCase.verifyWarning(@() logger.warning('Lost.'), 'rooftopMoM:logging:sinkFailed');
        end

        function testResetRestoresTheDefaults(testCase)
            logger = rooftopMoM.logging.Logger.instance();
            logger.configure(rooftopMoM.logging.LogLevel.TRACE, {});

            logger.reset();

            testCase.verifyEqual(logger.level(), rooftopMoM.logging.LogLevel.WARNING);
        end

        function testDisabledCallIsCheap(testCase)
            logger = rooftopMoM.logging.Logger.instance();
            callCount = 1e5;

            start = tic;
            for n = 1:callCount
                logger.trace('Iteration %d.', n);
            end
            secondsPerCall = toc(start) / callCount;

            % A solver iteration takes milliseconds; a disabled call must stay below 20 microseconds.
            testCase.verifyLessThan(secondsPerCall, 20e-6);
        end
    end

    methods (Access = private)
        function sink = captureLines_(~, level)
            sink = rooftopMoM.logging.sinks.MemorySink();
            rooftopMoM.logging.Logger.instance().configure(level, {sink});
        end
    end
end

function warnFromLocalFunction_()
    rooftopMoM.logging.Logger.instance().warning('From a local function.');
end
