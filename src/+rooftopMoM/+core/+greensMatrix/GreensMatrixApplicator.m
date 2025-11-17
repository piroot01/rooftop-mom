classdef GreensMatrixApplicator < handle
    properties (Access = protected)
        context_;
        matrixData_ rooftopMoM.core.greensMatrix.MatrixData;
    end

    methods
        function obj = GreensMatrixApplicator(context, matrixData)
            obj.context_ = context;
            obj.matrixData_ = matrixData;
        end
    end

    methods (Abstract)
        generate(obj);

        output = apply(obj, input);

        varargout = inputSize(obj);

        memory = memory(obj);
    end
end
