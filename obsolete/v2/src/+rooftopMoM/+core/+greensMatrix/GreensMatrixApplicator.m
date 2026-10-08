classdef GreensMatrixApplicator < handle
% GreensMatrixApplicator: Base of the matrix-free Green's matrix applicators that work with spectral matrix data.

    properties (Access = protected)
        context_;
        spectralMatrixData_ rooftopMoM.core.greensMatrix.SpectralMatrixData;
        rooftopContent_ rooftopMoM.basisFunctions.RooftopContent;
    end

    methods
        function obj = GreensMatrixApplicator(context, spectralMatrixData, rooftopContent)
            obj.context_ = context;
            obj.spectralMatrixData_ = spectralMatrixData;
            obj.rooftopContent_ = rooftopContent;
        end
    end

    methods
        function varargout = inputSize(obj)
            if nargout <= 1
                varargout{1} = obj.context_.rooftops().count();
            elseif nargout == 2
                varargout{1} = obj.context_.rooftops().count();
                varargout{2} = 1;
            end
        end

        function memory = memory(obj)
        % memory: Returns the memory occupied by the stored spectral matrix data.
        %
        %   Output:
        %       memory: rooftopMoM.Memory of the spectral blocks.

            memory = obj.spectralMatrixData_.memory();
        end

        function rooftopContent = rooftopContent(obj)
            rooftopContent = obj.rooftopContent_;
        end
    end

    methods (Abstract)
        generate(obj);

        output = apply(obj, input);

        % workingMemory: Returns the peak memory of the arrays alive during one application (output included).
        memory = workingMemory(obj);
    end
end
