classdef FFTAcceleratedApplicator < rooftopMoM.core.greensMatrix.GreensMatrixApplicator
    methods
        function obj = FFTAcceleratedApplicator(context, spectralMatrixData, rooftopContent)
            obj@rooftopMoM.core.greensMatrix.GreensMatrixApplicator(context, spectralMatrixData, rooftopContent);
        end
    end

    methods
        function generate(~)
            % NOTE: Nothing has to be prepared
        end

        function memory = workingMemory(obj)
        % workingMemory: Returns the peak memory of the arrays alive during one application (output included).
        %
        %   Output:
        %       memory: rooftopMoM.Memory of the peak working arrays.

            [blockSizeX, blockSizeY] = obj.spectralMatrixData_.blockSize();
            complexBytes = rooftopMoM.DataTypeMemory.memory(rooftopMoM.DataTypeMemory.ComplexDouble).bytes();
            blockBytes = blockSizeX * blockSizeY * complexBytes;

            switch obj.rooftopContent_
                case rooftopMoM.basisFunctions.RooftopContent.XONLY
                    % Output spectrum and ifft2 of the n-by-2 block, which needs one and a half blocks.
                    memory = rooftopMoM.Memory(2.5 * blockBytes, 'B');
                case rooftopMoM.basisFunctions.RooftopContent.YONLY
                    % Output spectrum and ifft2 of the 2-by-n block, which needs two blocks.
                    memory = rooftopMoM.Memory(3 * blockBytes, 'B');
                case rooftopMoM.basisFunctions.RooftopContent.BOTH
                    % Both output spectra, the padded spectrum of the y-directed input and its slice.
                    primitive = obj.context_.rooftops().primitive();
                    sliceBytes = double(primitive.countX()) * double(primitive.countY() - 1) * complexBytes;
                    memory = rooftopMoM.Memory(3 * blockBytes + sliceBytes, 'B');
                otherwise
                    error('Unknown rooftop content: %s', char(obj.rooftopContent_));
            end
        end

        function output = apply(obj, input)
            basisSizeX = obj.context_.rooftops().countDir('x');
            basisSizeY = obj.context_.rooftops().countDir('y');

            switch obj.rooftopContent_
                case rooftopMoM.basisFunctions.RooftopContent.XONLY
                    % X-only basis functions: inputs are already in the
                    % expected vector layout (no explicit reshape here).
                    % Zero-extend and FFT.
                    [blockSizeX, blockSizeY] = obj.spectralMatrixData_.blockSize();
                    inputSpectr = fft2(input, blockSizeX, blockSizeY);

                    % Spectral multiplication: G_xx * FFT(input).
                    outputSpectr = obj.spectralMatrixData_.blockXX() .* inputSpectr;
                    clear inputSpectr;

                    % Inverse FFT back to spatial domain.
                    output = ifft2(outputSpectr);
                    clear outputSpectr;

                    % Truncate the useful part and return as a column (truncate returns an array with shape [sizeX, 1]).
                    output = output(1:basisSizeX, 1);
                case rooftopMoM.basisFunctions.RooftopContent.YONLY 
                    % Y-only basis functions: the physical layout requires transposition before padding/FFT.
                    input = input.';

                    % Zero-extend and FFT (transpose handled above).
                    [blockSizeX, blockSizeY] = obj.spectralMatrixData_.blockSize();
                    inputSpectr = fft2(input, blockSizeX, blockSizeY);

                    % Spectral multiplication.
                    outputSpectr = obj.spectralMatrixData_.blockYY() .* inputSpectr;
                    clear inputSpectr;

                    % Inverse FFT.
                    output = ifft2(outputSpectr);
                    clear outputSpectr;

                    % Truncate the useful region and return vectorized.
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

                    outputXSpectr = obj.spectralMatrixData_.blockXX() .* inputXSpectr;
                    outputYSpectr = obj.spectralMatrixData_.blockYX() .* inputXSpectr;
                    clear inputXSpectr;

                    % Secondly we process the y direction of the input
                    inputYReshaped = reshape(input((basisSizeX + 1):end), countX, countY - 1);

                    inputYSpectr = fft2(inputYReshaped, blockSizeX, blockSizeY);
                    clear inputYReshaped;

                    outputXSpectr = outputXSpectr + obj.spectralMatrixData_.blockXY() .* inputYSpectr;
                    outputYSpectr = outputYSpectr + obj.spectralMatrixData_.blockYY() .* inputYSpectr;
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
    end
end
