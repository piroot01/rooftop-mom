classdef (Abstract) Base < matlab.unittest.TestCase
% Base: Common fixture of the unit tests. It puts src on the path for the duration of the test class, so a test class
%   can also be run on its own with only the folder tests on the path.

    methods (TestClassSetup)
        function addSourcePath(testCase)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'src')));
        end
    end
end
