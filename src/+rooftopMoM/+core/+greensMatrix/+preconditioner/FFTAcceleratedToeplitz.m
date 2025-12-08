classdef FFTAcceleratedToeplitz < rooftopMoM.core.greensMatrix.preconditioner.PreconditionerInterface
    properties (Access = private)
        spectralMatrixData_ rooftopMoM.core.greensMatrix.SpectralMatrixData;

        determinant_;
    end

    methods
        function obj = FFTAcceleratedToeplitz(context, spectralMatrixData)
            obj@rooftopMoM.core.greensMatrix.preconditioner.PreconditionerInterface(context);

            obj.spectralMatrixData_ = spectralMatrixData;
        end
    end

    methods
        function generate(obj)
            obj.determinant_ = rooftopMoM.core.greensMatrix.preconditioner.utils.generateSpectralGreensMatrixDeterminant( ...
                obj.spectralMatrixData_.gxxf(), ...
                obj.spectralMatrixData_.gxyf(), ...
                obj.spectralMatrixData_.gyxf(), ...
                obj.spectralMatrixData_.gyyf(), ...
                obj.spectralMatrixData_.rooftopContent() ...
            );
        end

        function output = apply(obj, input)
            sizeX = obj.context_.rooftops().countDir('x');
            sizeY = obj.context_.rooftops().countDir('y');

            switch obj.spectralMatrixData_.rooftopContent()
                case rooftopMoM.basisFunctions.RooftopContent.XONLY
                    [circSizeX, circSizeY] = size(obj.spectralMatrixData_.gxxf());
                    inputSpectr = fft2(input, circSizeX, circSizeY);

                    outputSpectr = inputSpectr ./ obj.spectralMatrixData_.gxxf();
                    clear inputSpectr;

                    output = ifft2(outputSpectr);
                    clear outputSpectr;

                    output = output(1:sizeX, 1);
                case rooftopMoM.basisFunctions.RooftopContent.YONLY
                    input = input.';

                    [circSizeX, circSizeY] = size(obj.spectralMatrixData_.gyyf());
                    inputSpectr = fft2(input, circSizeX, circSizeY);

                    outputSpectr = inputSpectr ./ obj.spectralMatrixData_.gyyf();
                    clear inputSpectr;

                    % Inverse FFT.
                    output = ifft2(outputSpectr);
                    clear outputSpectr;

                    output = output(1, 1:sizeY).';
                case rooftopMoM.basisFunctions.RooftopContent.BOTH
                    primitive = obj.context_.rooftops().primitive();
                    countX = primitive.countX();
                    countY = primitive.countY();

                    % First we process the x direction of the input
                    inputXReshaped = reshape(input(1:sizeX), countX - 1, countY);

                    [circSizeX, circSizeY] = size(obj.spectralMatrixData_.gxxf());
                    inputXSpectr = fft2(inputXReshaped, circSizeX, circSizeY);
                    clear inputXReshaped;

                    outputXSpectr = obj.spectralMatrixData_.gyyf() .* inputXSpectr;
                    outputYSpectr = (-1) * obj.spectralMatrixData_.gyxf() .* inputXSpectr;
                    clear inputXSpectr;

                    % Secondly we process the y direction of the input
                    inputYReshaped = reshape(input((sizeX + 1):end), countX, countY - 1);

                    inputYSpectr = fft2(inputYReshaped, circSizeX, circSizeY);
                    clear inputYReshaped;

                    outputXSpectr = outputXSpectr - obj.spectralMatrixData_.gxyf() .* inputYSpectr;
                    outputYSpectr = outputYSpectr + obj.spectralMatrixData_.gxxf() .* inputYSpectr;
                    clear inputYSpectr;

                    outputXSpectr = outputXSpectr ./ obj.determinant_;
                    outputYSpectr = outputYSpectr ./ obj.determinant_;

                    % And apply the inverse FFTs sequentially
                    outputX = ifft2(outputXSpectr);
                    clear outputXSpectr; % spatial domain

                    outputXTruncated = outputX(1:(countX - 1), 1:countY);
                    clear outputX;

                    outputY = ifft2(outputYSpectr);
                    clear outputYSpectr;

                    outputYTruncated = outputY(1:countX, 1:(countY - 1));
                    clear outputY;

                    % final concatenated output
                    output = [outputXTruncated(:); outputYTruncated(:)];
                otherwise
                    error('Unknown rooftop content: %s', rooftopContent);
            end
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
