classdef FFTAcceleratedApplicator < rooftopMoM.core.greensMatrix.GreensMatrixApplicator
    properties (Access = private)
        gxxf_;
        gxyf_;
        gyxf_;
        gyyf_;

        determinant_;
    end

    methods
        function obj = FFTAcceleratedApplicator(context, matrixData)
            obj@rooftopMoM.core.greensMatrix.GreensMatrixApplicator(context, matrixData);
        end
    end

    methods
        function generate(obj)
            [magneticXXCirc, magneticYYCirc, electricCirc] = rooftopMoM.core.greensMatrix.utils.generateGreensMatrixCirculants(obj.matrixData_);

            primitive = obj.context_.rooftops().primitive();
            [obj.gxxf_, obj.gxyf_, obj.gyxf_, obj.gyyf_, obj.determinant_] = rooftopMoM.core.greensMatrix.utils.generateSpectralGreensMatrix( ...
                magneticXXCirc, magneticYYCirc, electricCirc, ...
                obj.matrixData_.rooftopContent(), ...
                primitive.sizeX(), ...
                primitive.sizeY() ...
            );
        end

        function output = apply(obj, input)
            sizeX = obj.context_.rooftops().countDir('x');
            sizeY = obj.context_.rooftops().countDir('y');

            switch obj.matrixData_.rooftopContent()
                case rooftopMoM.basisFunctions.RooftopContent.XONLY
                    % X-only basis functions: inputs are already in the
                    % expected vector layout (no explicit reshape here).
                    % Zero-extend and FFT.
                    [circSizeX, circSizeY] = size(obj.gxxf_);
                    inputSpectr = fft2(input, circSizeX, circSizeY);

                    % Spectral multiplication: G_xx * FFT(input).
                    outputSpectr = obj.gxxf_ .* inputSpectr;
                    clear inputSpectr;

                    % Inverse FFT back to spatial domain.
                    output = ifft2(outputSpectr);
                    clear outputSpectr;

                    % Truncate the useful part and return as a column (truncate returns an array with shape [sizeX, 1]).
                    output = output(1:sizeX, 1);
                case rooftopMoM.basisFunctions.RooftopContent.YONLY 
                    % Y-only basis functions: the physical layout requires transposition before padding/FFT.
                    input = input.';

                    % Zero-extend and FFT (transpose handled above).
                    [circSizeX, circSizeY] = size(obj.gyyf_);
                    inputSpectr = fft2(input, circSizeX, circSizeY);

                    % Spectral multiplication.
                    outputSpectr = obj.gyyf_ .* inputSpectr;
                    clear inputSpectr;

                    % Inverse FFT.
                    output = ifft2(outputSpectr);
                    clear outputSpectr;

                    % Truncate the useful region and return vectorized.
                    output = output(1, 1:sizeY).';
                case rooftopMoM.basisFunctions.RooftopContent.BOTH
                    primitive = obj.context_.rooftops().primitive();
                    countX = primitive.countX();
                    countY = primitive.countY();

                    % First we process the x direction of the input
                    inputXReshaped = reshape(input(1:sizeX), countX - 1, countY);

                    [circSizeX, circSizeY] = size(obj.gxxf_);
                    inputXSpectr = fft2(inputXReshaped, circSizeX, circSizeY);
                    clear inputXReshaped;

                    outputXSpectr = obj.gxxf_ .* inputXSpectr;
                    outputYSpectr = obj.gyxf_ .* inputXSpectr;
                    clear inputXSpectr;

                    % Secondly we process the y direction of the input
                    inputYReshaped = reshape(input((sizeX + 1):end), countX, countY - 1);

                    inputYSpectr = fft2(inputYReshaped, circSizeX, circSizeY);
                    clear inputYReshaped;

                    outputXSpectr = outputXSpectr + obj.gxyf_ .* inputYSpectr;
                    outputYSpectr = outputYSpectr + obj.gyyf_ .* inputYSpectr;
                    clear inputYSpectr;

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

        function output = applyPreconditioner(obj, input)
            sizeX = obj.context_.rooftops().countDir('x');
            sizeY = obj.context_.rooftops().countDir('y');

            switch obj.matrixData_.rooftopContent()
                case rooftopMoM.basisFunctions.RooftopContent.XONLY
                    [circSizeX, circSizeY] = size(obj.gxxf_);
                    inputSpectr = fft2(input, circSizeX, circSizeY);

                    outputSpectr = inputSpectr ./ obj.gxxf_;
                    clear inputSpectr;

                    output = ifft2(outputSpectr);
                    clear outputSpectr;

                    output = output(1:sizeX, 1);
                case rooftopMoM.basisFunctions.RooftopContent.YONLY
                    input = input.';

                    [circSizeX, circSizeY] = size(obj.gyyf_);
                    inputSpectr = fft2(input, circSizeX, circSizeY);

                    outputSpectr = inputSpectr ./ obj.gyyf_;
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

                    [circSizeX, circSizeY] = size(obj.gxxf_);
                    inputXSpectr = fft2(inputXReshaped, circSizeX, circSizeY);
                    clear inputXReshaped;

                    outputXSpectr = obj.gyyf_ .* inputXSpectr;
                    outputYSpectr = (-1) * obj.gyxf_ .* inputXSpectr;
                    clear inputXSpectr;

                    % Secondly we process the y direction of the input
                    inputYReshaped = reshape(input((sizeX + 1):end), countX, countY - 1);

                    inputYSpectr = fft2(inputYReshaped, circSizeX, circSizeY);
                    clear inputYReshaped;

                    outputXSpectr = outputXSpectr - obj.gxyf_ .* inputYSpectr;
                    outputYSpectr = outputYSpectr + obj.gxxf_ .* inputYSpectr;
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
