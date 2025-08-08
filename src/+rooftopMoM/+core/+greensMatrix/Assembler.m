classdef Assembler < handle
    methods (Abstract)
        G = generateFull(obj, context);

        G = generateSpecified(obj, context, chunks);
    end
end
