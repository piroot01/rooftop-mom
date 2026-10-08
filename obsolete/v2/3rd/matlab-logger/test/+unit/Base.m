classdef (Abstract) Base < matlab.unittest.TestCase
    methods (TestMethodSetup)
        function addSrcToPath(~)
            % Find the project root by walking up until we find "src" folder
            p = fileparts(mfilename('fullpath'));

            while ~isempty(p) && ~exist(fullfile(p,'src'),'dir')
                p = fileparts(p);
            end

            if isempty(p)
                error('Could not find project root containing "src". Run tests from project tree or adjust detection.');
            end

            srcDir = fullfile(p, 'src');
            if ~contains(path, srcDir)
                addpath(srcDir);
            end

            testDir = fullfile(p, 'test');
            if ~contains(path, testDir)
                addpath(testDir);
            end

            % Reset singletons/persistent state.
            clear classes; %#ok<CLCLS>
        end
    end
end
