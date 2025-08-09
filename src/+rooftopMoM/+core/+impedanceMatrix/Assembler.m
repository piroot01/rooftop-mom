classdef Assembler < handle
    methods (Abstract)
        output = generateFull(obj, context);

        output = generateSpecified(obj, context, chunks);
    end
end
