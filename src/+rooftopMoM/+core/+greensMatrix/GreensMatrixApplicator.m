classdef GreensMatrixApplicator < handle
    properties (Access = protected)
        context_;
        spectralMatrixData_ rooftopMoM.core.greensMatrix.SpectralMatrixData;
    end

    methods
        function obj = GreensMatrixApplicator(context, spectralMatrixData)
            obj.context_ = context;
            obj.spectralMatrixData_ = spectralMatrixData;
        end
    end

    methods (Abstract)
        generate(obj);

        output = apply(obj, input);

        varargout = inputSize(obj);

        memory = memory(obj);
    end
end
