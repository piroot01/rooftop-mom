classdef PreconditionerInterface < handle
    properties (Access = protected)
        context_;
    end

    methods
        function obj = PreconditionerInterface(context)
            obj.context_ = context;
        end
    end

    methods (Abstract)
        output = apply(obj, input);

        varargout = inputSize(obj);

        memory = memory(obj);
    end
end
