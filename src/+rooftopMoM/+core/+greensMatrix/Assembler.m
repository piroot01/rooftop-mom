classdef Assembler < handle
    methods (Abstract)
        matrixData = generate(obj, context);

        G = generateFull(obj, context);

        G = generateSpecified(obj, context, chunks);
    end
end
