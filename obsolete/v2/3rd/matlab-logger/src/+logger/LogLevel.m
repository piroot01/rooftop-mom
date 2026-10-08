classdef LogLevel
    properties (Constant)
        TRACE = 0;
        DEBUG = 10;
        INFO = 20;
        WARN = 30;
        ERROR = 40;
        FATAL = 50;
        OFF = 1e9;
    end

    methods (Static)
        function level = fromName(name)
            if isempty(name)
                error('logger.LogLevel.fromName: empty name');
            end

            switch lower(string(name))
                case "trace", level = logger.LogLevel.TRACE;
                case "debug", level = logger.LogLevel.DEBUG;
                case "info", level = logger.LogLevel.INFO;
                case "warn", level = logger.LogLevel.WARN;
                case "error", level = logger.LogLevel.ERROR;
                case "fatal", level = logger.LogLevel.FATAL;
                case "off", level = logger.LogLevel.OFF;
                otherwise
                    error('Unknown log level name: %s', name);
            end
        end

        function name = name(level)
            if level == logger.LogLevel.TRACE
                name = "trace";
                return;
            end

            if level == logger.LogLevel.DEBUG
                name = "debug";
                return;
            end

            if level == logger.LogLevel.INFO
                name = "info";
                return;
            end

            if level == logger.LogLevel.WARN
                name = "warn";
                return;
            end

            if level == logger.LogLevel.ERROR
                name = "error";
                return;
            end

            if level == logger.LogLevel.FATAL
                name = "fatal";
                return;
            end

            if level == logger.LogLevel.OFF
                name = "off";
                return;
            end

            name = "level(" + string(level) + ")";
        end
    end
end
