classdef FFTAcceleratedToeplitz < rooftopMoM.core.greensMatrix.preconditioner.PreconditionerInterface
    properties (Access = private)
        rooftopContent_ rooftopMoM.basisFunctions.RooftopContent;

        determinant_;
    end

    methods
        function obj = FFTAcceleratedToeplitz(context, spectralMatrixData, rooftopContent)
            obj@rooftopMoM.core.greensMatrix.preconditioner.PreconditionerInterface(context, spectralMatrixData);

            obj.rooftopContent_ = rooftopContent;
        end
    end

    methods
        function generate(obj)
            determinant = obj.spectralMatrixData_.generateDeterminant(obj.rooftopContent_);

            if any(abs(determinant) < eps, 'all')
                error('Some values of the determinant are smaller than eps!');
            else
                obj.determinant_ = determinant;
            end
        end

        function output = apply(obj, input)
            basisSizeX = obj.context_.rooftops().countDir('x');
            basisSizeY = obj.context_.rooftops().countDir('y');

            switch obj.rooftopContent_
                case rooftopMoM.basisFunctions.RooftopContent.XONLY
                    [blockSizeX, blockSizeY] = obj.spectralMatrixData_.blockSize();
                    inputSpectr = fft2(input, blockSizeX, blockSizeY);

                    outputSpectr = inputSpectr ./ obj.spectralMatrixData_.blockXX();
                    clear inputSpectr;

                    output = ifft2(outputSpectr);
                    clear outputSpectr;

                    output = output(1:basisSizeX, 1);
                case rooftopMoM.basisFunctions.RooftopContent.YONLY
                    input = input.';

                    [blockSizeX, blockSizeY] = obj.spectralMatrixData_.blockSize();
                    inputSpectr = fft2(input, blockSizeX, blockSizeY);

                    outputSpectr = inputSpectr ./ obj.spectralMatrixData_.blockYY();
                    clear inputSpectr;

                    % Inverse FFT.
                    output = ifft2(outputSpectr);
                    clear outputSpectr;

                    output = output(1, 1:basisSizeY).';
                case rooftopMoM.basisFunctions.RooftopContent.BOTH
                    primitive = obj.context_.rooftops().primitive();
                    countX = primitive.countX();
                    countY = primitive.countY();

                    % First we process the x direction of the input
                    inputXReshaped = reshape(input(1:basisSizeX), countX - 1, countY);

                    [blockSizeX, blockSizeY] = obj.spectralMatrixData_.blockSize();
                    inputXSpectr = fft2(inputXReshaped, blockSizeX, blockSizeY);
                    clear inputXReshaped;

                    outputXSpectr = obj.spectralMatrixData_.blockYY() .* inputXSpectr;
                    outputYSpectr = (-1) * obj.spectralMatrixData_.blockYX() .* inputXSpectr;
                    clear inputXSpectr;

                    % Secondly we process the y direction of the input
                    inputYReshaped = reshape(input((basisSizeX + 1):end), countX, countY - 1);

                    inputYSpectr = fft2(inputYReshaped, blockSizeX, blockSizeY);
                    clear inputYReshaped;

                    outputXSpectr = outputXSpectr - obj.spectralMatrixData_.blockXY() .* inputYSpectr;
                    outputYSpectr = outputYSpectr + obj.spectralMatrixData_.blockXX() .* inputYSpectr;
                    clear inputYSpectr;

                    determinant = obj.spectralMatrixData_.generateDeterminant(obj.rooftopContent_);

                    outputXSpectr = outputXSpectr ./ determinant;
                    outputYSpectr = outputYSpectr ./ determinant;

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
