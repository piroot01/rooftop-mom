classdef LoggerTest < unit.Base
    methods(Test)
        function testLoggingToCaptureSink(tc)
            sink = unit.helpers.logger.sink.TestCaptureSink();
            logger.Logger.configure('enabled', true, 'minLevel', logger.LogLevel.DEBUG, 'sinks', { sink } );

            logger.Logger.log(logger.LogLevel.DEBUG, 'hello %d', 1);

            tc.verifyNotEmpty(sink.messages);
            msg = sink.messages{1};

            % Must contain the formatted message and the level name.
            tc.verifyThat(msg, matlab.unittest.constraints.ContainsSubstring('hello 1'));
            tc.verifyThat(lower(msg), matlab.unittest.constraints.ContainsSubstring('[debug]'));
        end

        function testDisabledLoggerIsNoOp(tc)
            sink = unit.helpers.logger.sink.TestCaptureSink();
            logger.Logger.configure('enabled', false, 'minLevel', logger.LogLevel.TRACE, 'sinks', { sink } );

            logger.Logger.info('won''t be stored');
            tc.verifyEmpty(sink.messages, 'When logger is disabled, sinks should not receive messages.');
        end

        function testMinLevelFiltering(tc)
            sink = unit.helpers.logger.sink.TestCaptureSink();
            logger.Logger.configure('enabled', true, 'minLevel', logger.LogLevel.INFO, 'sinks', { sink } );

            logger.Logger.log(logger.LogLevel.DEBUG, 'debug should not appear');
            logger.Logger.log(logger.LogLevel.INFO, 'info should appear');

            tc.verifyEqual(numel(sink.messages), 1);
            tc.verifyThat(sink.messages{1}, matlab.unittest.constraints.ContainsSubstring('info should appear'));
        end

        function testAddSinkAndBothReceive(tc)
            sink1 = unit.helpers.logger.sink.TestCaptureSink();
            sink2 = unit.helpers.logger.sink.TestCaptureSink();

            logger.Logger.configure('enabled', true, 'minLevel', logger.LogLevel.TRACE, 'sinks', { sink1 } );
            logger.Logger.addSink(sink2);

            logger.Logger.info('multi-sink test');
            tc.verifyEqual(numel(sink1.messages), 1);
            tc.verifyEqual(numel(sink2.messages), 1);
        end

        function testLoggerDoesNotThrowWhenOneSinkFails(tc)
            good = unit.helpers.logger.sink.TestCaptureSink();
            bad = unit.helpers.logger.sink.FaultySink();

            logger.Logger.configure('enabled', true, 'minLevel', logger.LogLevel.TRACE, 'sinks', { bad, good } );

            logger.Logger.info('survive faulty sink');
            tc.verifyEqual(numel(good.messages), 1, 'Good sink should still receive message even if other sink fails');
        end
    end
end
