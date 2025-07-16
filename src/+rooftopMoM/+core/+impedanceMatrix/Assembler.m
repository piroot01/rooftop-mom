classdef Assembler < handle
    methods (Abstract)
        Z = generateFull(obj, context);

        Z = generateSpecified(obj, context, chunks);
    end
end
