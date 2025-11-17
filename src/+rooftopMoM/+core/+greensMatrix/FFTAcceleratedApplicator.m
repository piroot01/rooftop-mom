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
            primitive = obj.context_.rooftops().primitive();
            [gridSizeX, gridSizeY] = primitive.gridSize();

            [magneticXXCirc, magneticYYCirc, electricCirc] = rooftopMoM.core.greensMatrix.utils.generateGreensMatrixCirculants(obj.matrixData_, gridSizeX + 1, gridSizeY + 1);

            primitive = obj.context_.rooftops().primitive();
            [obj.gxxf_, obj.gxyf_, obj.gyxf_, obj.gyyf_, obj.determinant_] = rooftopMoM.core.greensMatrix.utils.generateSpectralGreensMatrix( ...
                magneticXXCirc, magneticYYCirc, electricCirc, ...
                obj.matrixData_.rooftopContent(), ...
                primitive.sizeX(), ...
                primitive.sizeY() ...
            );
        end

        function output = apply(obj, input)
            % all blocks of the greens matrix have same size
            [circSizeX, circSizeY] = size(obj.gxxf_);

            sizeX = obj.context_.rooftops().countDir('x');
            sizeY = obj.context_.rooftops().countDir('y');

            switch obj.matrixData_.rooftopContent()
                case rooftopMoM.basisFunctions.RooftopContent.XONLY
                    % X-only basis functions: inputs are already in the
                    % expected vector layout (no explicit reshape here).
                    % Zero-extend and FFT.
                    inputSpectr = fft2(input, circSizeX, circSizeY);

                    % Spectral multiplication: G_xx * FFT(input).
                    outputSpectr = obj.gxxf_ .* inputSpectr;

                    % Inverse FFT back to spatial domain.
                    output = ifft2(outputSpectr);

                    % Truncate the useful part and return as a column (truncate returns an array with shape [sizeX, 1]).
                    output = output(1:sizeX, 1);
                case rooftopMoM.basisFunctions.RooftopContent.YONLY 
                    % Y-only basis functions: the physical layout requires transposition before padding/FFT.
                    input = input.';

                    % Zero-extend and FFT (transpose handled above).
                    inputSpectr = fft2(input, circSizeX, circSizeY);

                    % Spectral multiplication.
                    outputSpectr = obj.gyyf_ .* inputSpectr;

                    % Inverse FFT.
                    output = ifft2(outputSpectr);

                    % Truncate the useful region and return vectorized.
                    outputTruncated = output(1, 1:sizeY);

                    output = outputTruncated(:);
                case rooftopMoM.basisFunctions.RooftopContent.BOTH
                    % Old version
                    %{
                    % Both X and Y basis functions are present. The vector is split into two segments: 
                    % first the X-directed rooftop coefficients, then the Y-directed ones.
                    primitive = obj.context_.rooftops().primitive();
                    countX = primitive.countX();
                    countY = primitive.countY();

                    % Split vector into two parts corresponding to the two distinct grid samplings for X and Y roofops
                    inputX = input(1:sizeX);
                    inputY = input((sizeX + 1):end);

                    % Reshape into 2D arrays matching their grid layout.
                    inputXReshaped = reshape(inputX, countX - 1, countY);
                    inputYReshaped = reshape(inputY, countX, countY - 1);
                    clear inputX inputY;

                    % FFT each component after padding.
                    inputXSpectr = fft2(inputXReshaped, circSizeX, circSizeY);
                    inputYSpectr = fft2(inputYReshaped, circSizeX, circSizeY);
                    clear inputXReshaped inputYReshaped;

                    % Apply full 2x2 spectral Green's tensor.
                    outputX = ifft2(obj.gxxf_ .* inputXSpectr + obj.gxyf_ .* inputYSpectr);
                    outputY = ifft2(obj.gyxf_ .* inputXSpectr + obj.gyyf_ .* inputYSpectr);
                    clear inputXSpectr inputYSpectr;

                    % Truncate to useful region (removing padding) and re-vectorize both parts to produce final column.
                    outputXTruncated = outputX(1:(countX - 1), 1:countY);
                    outputYTruncated = outputY(1:countX, 1:(countY - 1));
                    clear outputX outputY;

                    output = [outputXTruncated(:); outputYTruncated(:)];
                    %}

                    %
                    % Memory optimised version
                    %

                    primitive = obj.context_.rooftops().primitive();
                    countX = primitive.countX();
                    countY = primitive.countY();

                    % First we process the x direction of the input
                    inputXReshaped = reshape(input(1:sizeX), countX - 1, countY);

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
