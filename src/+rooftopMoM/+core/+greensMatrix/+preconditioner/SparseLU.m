classdef SparseLU < rooftopMoM.core.greensMatrix.preconditioner.PreconditionerInterface
    properties (Access = private)
        count_;
        middleMode_;
        bandwidth_;

        L_;
        U_;
    end

    methods
        function obj = SparseLU(context, bandwidth)
            obj@rooftopMoM.core.greensMatrix.preconditioner.PreconditionerInterface(context);

            count = context.rooftops().count();

            middleMode = floor(count / 2);

            if (middleMode - bandwidth < 1) || (middleMode + bandwidth > count)
                error('Unknown count N = %d is too small for bandwidth = %d. Increase countX or reduce bandwidth.', count, bandwidth);
            end

            obj.count_ = count;
            obj.middleMode_ = middleMode;
            obj.bandwidth_ = bandwidth;
        end
    end

    methods
        function probeVector = generateProbeVector(obj)
            probeVector = zeros(obj.count_, 1);
            probeVector(obj.middleMode_) = 1;
        end

        function generate(obj, probeResponse)
            extractionIndices = (obj.middleMode_ - obj.bandwidth_):(obj.middleMode_ + obj.bandwidth_);
            extractedValues = probeResponse(extractionIndices);

            diagonalsIndices = obj.bandwidth_:(-1):((-1) * obj.bandwidth_);
            diagonals = repmat(extractedValues(:).', obj.count_, 1);

            preconditioner = spdiags(diagonals, diagonalsIndices, obj.count_, obj.count_);
            [obj.L_, obj.U_] = lu(preconditioner);
        end

        function output = apply(obj, input)
            output = obj.U_ \ (obj.L_ \ input);
        end

        function varargout = inputSize(obj)
            if nargout <= 1
                varargout{1} = obj.context_.rooftops().count();
            elseif nargout == 2
                varargout{1} = obj.context_.rooftops().count();
                varargout{2} = 1;
            end
        end

        function memory = memory(obj)
            % TODO: implement
        end
    end
end
