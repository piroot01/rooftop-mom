classdef FFTAcceleratedCirculant < rooftopMoM.core.greensMatrix.preconditioner.PreconditionerInterface
    properties (Access = private)
        spectralMatrixData_;
    end

    methods
        function obj = FFTAcceleratedCirculant(context, spectralMatrixData)
            obj@rooftopMoM.core.greensMatrix.preconditioner.PreconditionerInterface(context);

            obj.spectralMatrixData_ = spectralMatrixData;
        end
    end

    methods
        function output = apply(obj, input)
            sizeX = obj.context_.rooftops().countDir('x');
            sizeY = obj.context_.rooftops().countDir('y');

            switch obj.spectralMatrixData_.rooftopContent()
                case rooftopMoM.basisFunctions.RooftopContent.XONLY
                    [fftSizeX, ~] = size(obj.spectralMatrixData_.pxxf());
                    inputSpectr = fft(input, fftSizeX, 1);

                    outputSpectr = inputSpectr ./ obj.spectralMatrixData_.pxxf();
                    clear inputSpectr;

                    output = ifft(outputSpectr, [], 1);

                    output = output(1:sizeX);
                case rooftopMoM.basisFunctions.RooftopContent.YONLY
                    input = input.';

                    [~, fftSizeY] = size(obj.spectralMatrixData_.pyyf());
                    inputSpectr = fft(input, fftSizeY, 2);

                    outputSpectr = inputSpectr ./ obj.spectralMatrixData_.pyyf();
                    clear inputSpectr;

                    % Inverse FFT.
                    output = ifft(outputSpectr, [], 2);
                    clear outputSpectr;

                    output = output(1:sizeY).';
                case rooftopMoM.basisFunctions.RooftopContent.BOTH
                    primitive = obj.context_.rooftops().primitive();
                    countX = primitive.countX();
                    countY = primitive.countY();

                    % First we process the x direction of the input
                    inputXReshaped = reshape(input(1:sizeX), countX - 1, countY);

                    [fftSizeX, fftSizeY] = size(obj.spectralMatrixData_.pxxf());
                    inputXSpectr = fft2(inputXReshaped, fftSizeX, fftSizeY);
                    clear inputXReshaped;

                    outputXSpectr = obj.spectralMatrixData_.pyyf() .* inputXSpectr;
                    outputYSpectr = (-1) * obj.spectralMatrixData_.pyxf() .* inputXSpectr;
                    clear inputXSpectr;

                    % Secondly we process the y direction of the input
                    inputYReshaped = reshape(input((sizeX + 1):end), countX, countY - 1);

                    inputYSpectr = fft2(inputYReshaped, fftSizeX, fftSizeY);
                    clear inputYReshaped;

                    outputXSpectr = outputXSpectr - obj.spectralMatrixData_.pxyf() .* inputYSpectr;
                    outputYSpectr = outputYSpectr + obj.spectralMatrixData_.pxxf() .* inputYSpectr;
                    clear inputYSpectr;

                    outputXSpectr = outputXSpectr ./ obj.spectralMatrixData_.determinant();
                    outputYSpectr = outputYSpectr ./ obj.spectralMatrixData_.determinant();

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
