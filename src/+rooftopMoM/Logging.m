classdef Logging
% Logging: User control of the logging of the rooftopMoM package.
%   Write the messages up to DEBUG to the command window and append them to run.log:
%       rooftopMoM.Logging.configure('level', 'debug', 'file', 'run.log');
%   Restore the default (WARNING, command window):
%       rooftopMoM.Logging.reset();
%   rooftopMoM.logging.LogLevel describes what each level reports.

    methods (Static)
        function configure(options)
        % configure: Sets the level and the destinations of the log messages.
        %
        %   Optional name-value inputs:
        %       'level': 'off', 'error', 'warning', 'info', 'debug' or 'trace' (default 'warning').
        %       'file': Text file the lines are appended to (default none).
        %       'console': Logical flag, write the lines to the command window (default true).

            arguments
                options.level {mustBeTextScalar} = 'warning'
                options.file {mustBeTextScalar} = ''
                options.console (1, 1) logical = true
            end

            levelName = char(options.level);
            fileName = char(options.file);
            levelNames = {'off', 'error', 'warning', 'info', 'debug', 'trace'};

            if ~any(strcmpi(levelName, levelNames))
                error('rooftopMoM:logging:unknownLevel', 'Unknown log level ''%s''. Use one of: %s.', ...
                    levelName, strjoin(levelNames, ', '));
            end

            sinks = {};

            if options.console
                sinks{end + 1} = rooftopMoM.logging.sinks.ConsoleSink();
            end

            if ~isempty(fileName)
                rooftopMoM.Logging.checkWritable_(fileName);
                sinks{end + 1} = rooftopMoM.logging.sinks.FileSink(fileName);
            end

            rooftopMoM.logging.Logger.instance().configure(rooftopMoM.logging.LogLevel.(upper(levelName)), sinks);
        end

        function reset()
        % reset: Restores the default logging: level WARNING, written to the command window.

            rooftopMoM.logging.Logger.instance().reset();
        end
    end

    methods (Static, Access = private)
        function checkWritable_(fileName)
        % checkWritable_: Raises an error if the file cannot be opened for appending.
        %
        %   Input:
        %       fileName: Character vector of the file name.

            % Fail here once, instead of a warning on every logged line.
            [fileId, message] = fopen(fileName, 'a');

            if fileId < 0
                error('rooftopMoM:logging:fileNotWritable', 'Cannot open the log file %s: %s', fileName, message);
            end

            fclose(fileId);
        end
    end
end
