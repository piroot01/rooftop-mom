classdef FileSink < rooftopMoM.logging.sinks.Sink
% FileSink: Appends the log lines to a text file.
%   The file is opened and closed for each line. Every line is on the disk at once, and the file is never held
%   open, which on Windows would lock it.

    properties (Access = private)
        fileName_ char;
    end

    methods
        function obj = FileSink(fileName)
            obj.fileName_ = fileName;
        end
    end

    methods
        function write(obj, line)
        % write: Appends the line to the file.
        %
        %   Input:
        %       line: Character vector of the formatted line.

            [fileId, message] = fopen(obj.fileName_, 'a', 'n', 'UTF-8');

            if fileId < 0
                error('rooftopMoM:logging:fileNotWritable', 'Cannot open the log file %s: %s', ...
                    obj.fileName_, message);
            end

            fprintf(fileId, '%s\n', line);
            fclose(fileId);
        end

        function fileName = fileName(obj)
            fileName = obj.fileName_;
        end
    end
end
