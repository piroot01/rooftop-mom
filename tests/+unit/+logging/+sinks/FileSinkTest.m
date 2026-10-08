classdef FileSinkTest < unit.Base
% FileSinkTest: Tests of rooftopMoM.logging.sinks.FileSink.

    methods (Test)
        function testAppendsLines(testCase)
            fileName = testCase.temporaryFile_();
            sink = rooftopMoM.logging.sinks.FileSink(fileName);

            sink.write('first');
            sink.write('second');

            testCase.verifyEqual(cellstr(readlines(fileName, 'EmptyLineRule', 'skip')), {'first'; 'second'});
        end

        function testKeepsTheExistingContent(testCase)
            fileName = testCase.temporaryFile_();
            writelines('existing', fileName);

            rooftopMoM.logging.sinks.FileSink(fileName).write('appended');

            testCase.verifyEqual(cellstr(readlines(fileName, 'EmptyLineRule', 'skip')), {'existing'; 'appended'});
        end

        function testUnwritableFileRaisesError(testCase)
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture()).Folder;
            sink = rooftopMoM.logging.sinks.FileSink(folder);

            testCase.verifyError(@() sink.write('line'), 'rooftopMoM:logging:fileNotWritable');
        end
    end

    methods (Access = private)
        function fileName = temporaryFile_(testCase)
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture()).Folder;
            fileName = fullfile(folder, 'run.log');
        end
    end
end
