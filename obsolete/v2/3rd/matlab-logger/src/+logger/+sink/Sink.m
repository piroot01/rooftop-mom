classdef (Abstract) Sink < handle
    methods (Abstract)
        write(obj, formattedMessage);
    end
end
