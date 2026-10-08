classdef FileSinkTest < unit.Base
    properties (Access = private)
        tmpFileName_ = '';
    end

    methods(TestMethodSetup)
        function makeTmp(tc)
            tc.tmpFileName_ = [tempname() '.log'];

            % Ensure file does not exist yet.
            if exist(tc.tmpFileName_, 'file')
                delete(tc.tmpFileName_);
            end
        end
    end

    methods(TestMethodTeardown)
        function cleanup(tc)
            if exist(tc.tmpFileName_, 'file')
                delete(tc.tmpFileName_);
            end
        end
    end

    methods(Test)
        function testWriteCreatesAndAppends(tc)
            f = logger.sink.FileSink(tc.tmpFileName_);

            % Write two messages.
            f.write('[t1] message1\n');
            f.write('[t2] message2\n');

            % Close to flush.
            delete(f);

            data = fileread(tc.tmpFileName_);
            tc.verifyThat(data, matlab.unittest.constraints.ContainsSubstring('message1'));
            tc.verifyThat(data, matlab.unittest.constraints.ContainsSubstring('message2'));
        end

        function testRecoveryOpenOnWrite(tc)
            f = logger.sink.FileSink(tc.tmpFileName_);

            % Simulate close.
            f.close();

            % Still should reopen on write.
            f.write('[t3] afterclose\n');
            delete(f);

            data = fileread(tc.tmpFileName_);
            tc.verifyThat(data, matlab.unittest.constraints.ContainsSubstring('afterclose'));
        end
    end
end
