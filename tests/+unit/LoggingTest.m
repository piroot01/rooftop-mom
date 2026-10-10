classdef LoggingTest < unit.Base
% LoggingTest: Tests of rooftopMoM.Logging.

    methods (Test)
        function testWritesToTheFile(testCase)
            fileName = testCase.temporaryFile_();

            rooftopMoM.Logging.configure('level', 'debug', 'file', fileName, 'console', false);
            rooftopMoM.logging.Logger.instance().debug('Written.');

            lines = cellstr(readlines(fileName, 'EmptyLineRule', 'skip'));
            testCase.verifyNumElements(lines, 1);
            testCase.verifySubstring(lines{1}, '[debug]');
        end

        function testLevelIgnoresCase(testCase)
            rooftopMoM.Logging.configure('level', 'TRACE', 'console', false);
            testCase.verifyEqual(rooftopMoM.logging.Logger.instance().level(), rooftopMoM.logging.LogLevel.TRACE);
        end

        function testUnknownLevelRaisesError(testCase)
            testCase.verifyError(@() rooftopMoM.Logging.configure('level', 'verbose'), ...
                'rooftopMoM:logging:unknownLevel');
        end

        function testUnwritableFileRaisesError(testCase)
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture()).Folder;
            testCase.verifyError(@() rooftopMoM.Logging.configure('file', folder), ...
                'rooftopMoM:logging:fileNotWritable');
        end

        function testResetRestoresTheDefaults(testCase)
            rooftopMoM.Logging.configure('level', 'trace', 'console', false);
            rooftopMoM.Logging.reset();
            testCase.verifyEqual(rooftopMoM.logging.Logger.instance().level(), rooftopMoM.logging.LogLevel.WARNING);
        end
    end

    methods (Access = private)
        function fileName = temporaryFile_(testCase)
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture()).Folder;
            fileName = fullfile(folder, 'run.log');
        end
    end
end
