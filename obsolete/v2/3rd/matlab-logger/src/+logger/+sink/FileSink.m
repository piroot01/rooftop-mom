classdef FileSink < logger.sink.Sink
    properties (Access = private)
        fid_ = -1;
        filename_ (1, :) char = '';
    end

    methods
        function obj = FileSink(filename)
            if nargin == 0 || isempty(filename)
                error('logger.sink.FileSink.FileSink: empty file name');
            end

            obj.filename_ = filename;

            obj.open();
        end

        function delete(obj)
            obj.close();
        end

        function write(obj, formattedMessage)
            if obj.fid_ < 0
                obj.open();
            end

            if obj.fid_ >= 0
                fprintf(obj.fid_, '%s', formattedMessage);
            end
        end

        function open(obj)
            try
                [fid, msg] = fopen(obj.filename_, 'a');
                if fid < 0
                    warning('logger.sink.FileSink.open: cannot open file %s: %s', obj.filename, msg);
                else
                    obj.fid_ = fid;
                end
            catch ex
                warning('logger.sink.FileSink.open: failed: %s', ex.message());
            end
        end

        function close(obj)
            if obj.fid_ >= 0
                fclose(obj.fid_);
                obj.fid_ = -1;
            end
        end
    end
end
