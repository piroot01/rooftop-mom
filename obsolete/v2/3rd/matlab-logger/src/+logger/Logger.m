classdef Logger < handle
    properties (Access = private)
        sinks_ cell = {};
        minLevel_ double = logger.LogLevel.INFO;
        enabled_ logical = true;
        timeFormat_ char = 'yyyy-MM-dd HH:mm:ss.SSS';
        template_ char = '[%s] [%s] %s\n';
    end

    methods (Access = private)
        function obj = Logger()
            % The default: console sink.
            obj.sinks_ = { logger.sink.ConsoleSink() };
        end
    end

    methods (Static, Access = private)
        function instance = getInstance_()
            persistent INST

            if isempty(INST) || ~isvalid(INST)
                INST = logger.Logger();
            end

            instance = INST;
        end
    end

    methods (Access = private)
        function configure_(obj, varargin)
            % Example usage:
            %   configure('minLevel', logger.LogLevel.DEBUG, 'enabled', true, 'sinks', {logger.sink.ConsoleSink(), logger.sink.FileSink('a.log')})

            p = inputParser();
            addParameter(p, 'minLevel', obj.minLevel_);
            addParameter(p, 'enabled', obj.enabled_);
            addParameter(p, 'sinks', obj.sinks_);
            addParameter(p, 'timeFormat', obj.timeFormat_);
            addParameter(p, 'template', obj.template_);

            parse(p, varargin{:});

            obj.minLevel_ = p.Results.minLevel;
            obj.enabled_ = p.Results.enabled;
            obj.sinks_ = p.Results.sinks;
            obj.timeFormat_ = p.Results.timeFormat;
            obj.template_ = p.Results.template;
        end

        function addSink_(obj, sink)
            if isa(sink, 'logger.sink.Sink')
                obj.sinks_{end + 1} = sink;
            else
                error('logger.Logger.addSink: requires an instance of Sink');
            end
        end

        function setMinLevel_(obj, level)
            obj.minLevel_ = level;
        end

        function setEnabled_(obj, flag)
            obj.enabled_ = logical(flag);
        end

        function log_(obj, level, fmt, varargin)
            % Very fast early exit if logging disabled or level too low.
            if ~obj.enabled_ || level < obj.minLevel_
                return;
            end

            % Build timestamp only when needed.
            ts = char(datetime('now','Format', obj.timeFormat_));
            msg = sprintf(fmt, varargin{:});
            formatted = sprintf(obj.template_, ts, char(logger.LogLevel.name(level)), msg);

            % Send to sinks.
            for k=1:numel(obj.sinks_)
                try
                    obj.sinks_{k}.write(formatted);
                catch ex
                    % Do not throw from logger; optionally emit warning once.
                    warning('logger.Logger.log: sink write failed: %s', ex.message());
                end
            end
        end
    end

    methods (Static)
        function configure(varargin)
            inst = logger.Logger.getInstance_();
            inst.configure_(varargin{:});
        end

        function addSink(sink)
            inst = logger.Logger.getInstance_();
            inst.addSink_(sink);
        end

        function setMinLevel(level)
            inst = logger.Logger.getInstance_();
            inst.setMinLevel_(level);
        end

        function setEnabled(flag)
            inst = logger.Logger.getInstance_();
            inst.setEnabled_(flag);
        end

        function log(level, fmt, varargin)
            inst = logger.Logger.getInstance_();
            inst.log_(level, fmt, varargin{:});
        end

        function trace(fmt, varargin), logger.Logger.log(logger.LogLevel.TRACE, fmt, varargin{:}); end
        function debug(fmt, varargin), logger.Logger.log(logger.LogLevel.DEBUG, fmt, varargin{:}); end
        function info(fmt, varargin), logger.Logger.log(logger.LogLevel.INFO,  fmt, varargin{:}); end
        function warn(fmt, varargin), logger.Logger.log(logger.LogLevel.WARN,  fmt, varargin{:}); end
        function error(fmt, varargin), logger.Logger.log(logger.LogLevel.ERROR, fmt, varargin{:}); end
        function fatal(fmt, varargin), logger.Logger.log(logger.LogLevel.FATAL, fmt, varargin{:}); end
    end
end
