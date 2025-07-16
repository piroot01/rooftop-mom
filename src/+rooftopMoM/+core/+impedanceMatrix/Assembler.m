classdef Assembler < handle
    properties (Access = private)
        constants_ = rooftopMoM.core.Constants;
    end

    methods (Abstract)
        Z = generate(obj, context);
    end
end
