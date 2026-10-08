classdef PreconditionerInterface < handle
    properties (Access = protected)
        context_;
        spectralMatrixData_;
    end

    methods
        function obj = PreconditionerInterface(context, spectralMatrixData)
            obj.context_ = context;
            obj.spectralMatrixData_ = spectralMatrixData;
        end
    end

    methods (Abstract)
        output = apply(obj, input);

        varargout = inputSize(obj);

        memory = memory(obj);
    end
end
