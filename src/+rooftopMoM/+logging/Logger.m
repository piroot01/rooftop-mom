classdef Logger < handle
% Logger: Shared logger of the rooftopMoM package.
%   All modules write to one instance, Logger.instance(). A message is written only if its level is enabled; a
%   disabled call returns after a single comparison. A written line looks like
%       [2026-10-08 16:41:23.123] [INFO] [IterativeMethod.solve]: ConjugateGradient converged in 12 iterations.
%   The caller is read from the call stack, so a call site passes only the message:
%       logger = rooftopMoM.logging.Logger.instance();
%       logger.info('%s converged in %d iterations.', method, iterationCount);
%
%   NOTE: The user configures the logger through rooftopMoM.Logging; the modules only write to it. Workers of a
%       parallel pool have their own instance, so nothing is logged inside parfor.

    properties (Access = private)
        % Numeric value of the LogLevel, compared on every call.
        level_ double;

        % Cell of rooftopMoM.logging.sinks.Sink receiving every written line.
        sinks_ cell;

        % Wall clock at the last configuration and the tic taken at that moment.
        date_ char;
        secondsOfDay_ double;
        clock_ uint64;
    end

    methods (Access = private)
        function obj = Logger()
            obj.reset();
        end
    end

    methods
        function configure(obj, level, sinks)
        % configure: Sets the level and the sinks.
        %
        %   Inputs:
        %       level: rooftopMoM.logging.LogLevel of the most verbose messages that are written.
        %       sinks: Cell of rooftopMoM.logging.sinks.Sink receiving every written line.

            obj.level_ = double(level);
            obj.sinks_ = sinks;
            obj.readClock_();
        end

        function reset(obj)
        % reset: Restores the default configuration: level WARNING, written to the command window.

            obj.configure(rooftopMoM.logging.LogLevel.WARNING, {rooftopMoM.logging.sinks.ConsoleSink()});
        end

        function level = level(obj)
            level = rooftopMoM.logging.LogLevel(obj.level_);
        end

        function flag = isEnabled(obj, level)
        % isEnabled: Returns true if the messages of the given level are written. Use it to skip building a costly
        %   message, or read it once before a loop to guard the messages inside the loop.
        %
        %   Input:
        %       level: rooftopMoM.logging.LogLevel of the message.
        %
        %   Output:
        %       flag: Logical scalar.

            flag = double(level) <= obj.level_;
        end

        function error(obj, format, varargin)
        % error: Writes the message (sprintf format and values) at the level ERROR.

            if obj.level_ >= 1
                obj.write_('ERROR', format, varargin);
            end
        end

        function warning(obj, format, varargin)
        % warning: Writes the message (sprintf format and values) at the level WARNING.

            if obj.level_ >= 2
                obj.write_('WARNING', format, varargin);
            end
        end

        function info(obj, format, varargin)
        % info: Writes the message (sprintf format and values) at the level INFO.

            if obj.level_ >= 3
                obj.write_('INFO', format, varargin);
            end
        end

        function debug(obj, format, varargin)
        % debug: Writes the message (sprintf format and values) at the level DEBUG.

            if obj.level_ >= 4
                obj.write_('DEBUG', format, varargin);
            end
        end

        function trace(obj, format, varargin)
        % trace: Writes the message (sprintf format and values) at the level TRACE.

            if obj.level_ >= 5
                obj.write_('TRACE', format, varargin);
            end
        end
    end

    methods (Access = private)
        function write_(obj, levelName, format, values)
        % write_: Formats the line and passes it to every sink.
        %
        %   Inputs:
        %       levelName: Character vector of the level name.
        %       format: sprintf format of the message.
        %       values: Cell of the values of the format.

            % Skip the frames of write_ and of the level method; the first frame left is the caller.
            stack = dbstack(2, '-completenames');

            if isempty(stack)
                caller = 'command line';
            else
                caller = rooftopMoM.logging.Logger.callerName_(stack(1));
            end

            line = sprintf('[%s] [%s] [%s]: %s', obj.timestamp_(), levelName, caller, sprintf(format, values{:}));

            for n = 1:numel(obj.sinks_)
                % A failing sink must not stop the computation, so it only raises a warning.
                try
                    obj.sinks_{n}.write(line);
                catch exception
                    warning('rooftopMoM:logging:sinkFailed', 'A log sink failed: %s', exception.message);
                end
            end
        end

        function timestamp = timestamp_(obj)
        % timestamp_: Returns the current time as 'yyyy-MM-dd HH:mm:ss.SSS'.
        %
        %   Output:
        %       timestamp: Character vector of the time.

            % Formatting a datetime takes 0.4 ms, so the time is the clock read at configuration plus the tic time.
            milliseconds = floor(1000 * (obj.secondsOfDay_ + toc(obj.clock_)));

            % Read the clock again after midnight, so the date stays right.
            if milliseconds >= 86400000
                obj.readClock_();
                milliseconds = floor(1000 * obj.secondsOfDay_);
            end

            timestamp = sprintf('%s %02d:%02d:%02d.%03d', obj.date_, floor(milliseconds / 3600000), ...
                floor(mod(milliseconds, 3600000) / 60000), floor(mod(milliseconds, 60000) / 1000), ...
                mod(milliseconds, 1000));
        end

        function readClock_(obj)
        % readClock_: Reads the wall clock and starts the tic that the timestamps are measured from.

            moment = datetime('now');
            obj.clock_ = tic;
            obj.date_ = char(datetime(moment, 'Format', 'yyyy-MM-dd'));
            obj.secondsOfDay_ = seconds(timeofday(moment));
        end
    end

    methods (Static)
        function logger = instance()
        % instance: Returns the logger shared by the whole package, created with the default configuration on the
        %   first call.
        %
        %   Output:
        %       logger: rooftopMoM.logging.Logger.

            persistent shared

            if isempty(shared) || ~isvalid(shared)
                shared = rooftopMoM.logging.Logger();
            end

            logger = shared;
        end
    end

    methods (Static, Access = private)
        function caller = callerName_(frame)
        % callerName_: Returns the caller of a stack frame: ClassName.methodName for a method, the function name
        %   otherwise.
        %
        %   Input:
        %       frame: Element of the struct array returned by dbstack with the option '-completenames'.
        %
        %   Output:
        %       caller: Character vector of the caller name.

            [folder, fileName] = fileparts(frame.file);
            [~, folderName] = fileparts(folder);

            % dbstack names a method in an @ClassName folder by the method name only.
            if startsWith(folderName, '@')
                caller = [folderName(2:end), '.', frame.name];

            % dbstack names a local function by its own name only.
            elseif ~contains(frame.name, '.') && ~strcmp(frame.name, fileName)
                caller = [fileName, '>', frame.name];
            else
                caller = frame.name;
            end
        end
    end
end
