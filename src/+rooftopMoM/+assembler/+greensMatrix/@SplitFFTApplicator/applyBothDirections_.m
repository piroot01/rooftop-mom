function output = applyBothDirections_(obj, input)
% applyBothDirections_: Applies the Green's matrix of the x- and y-directed rooftops.
%
%   Input:
%       input: Nx1 vector of the rooftop coefficients, the x-directed rooftops first.
%
%   Output:
%       output: Nx1 vector of the product.
%
%   NOTE: The y direction is split first (2 FFTs per component) and the x direction second (4 FFTs per component along
%       the contiguous columns), which is the faster order. The arrays are cleared as soon as possible to keep the
%       operations in place.

    [rowCountX, columnCountX] = obj.rooftops_.gridSize('x');
    [rowCountY, columnCountY] = obj.rooftops_.gridSize('y');
    [fftSizeX, fftSizeY] = obj.spectra_(1).gridSize();
    countX = rowCountX * columnCountX;

    % Zero-extend both components to the subband grid.
    inputX = reshape(input(1:countX), rowCountX, columnCountX);
    inputX(fftSizeX, fftSizeY) = 0;

    inputY = reshape(input((countX + 1):end), rowCountY, columnCountY);
    inputY(fftSizeX, fftSizeY) = 0;

    for parityY = 0:1
        if parityY
            inputX = inputX .* obj.modulationY_;
            inputY = inputY .* obj.modulationY_;
        end

        spectrumX = fft(inputX, [], 2);
        spectrumY = fft(inputY, [], 2);

        if parityY
            clear inputX inputY;
        end

        for parityX = 0:1
            if parityX
                spectrumX = spectrumX .* obj.modulationX_;
                spectrumY = spectrumY .* obj.modulationX_;
            end

            subbandX = fft(spectrumX, [], 1);
            subbandY = fft(spectrumY, [], 1);

            if parityX
                clear spectrumX spectrumY;
            end

            blocks = obj.spectra_(1 + parityX + 2 * parityY);
            subbandOutputX = ifft(blocks.blockXX() .* subbandX + blocks.blockXY() .* subbandY, [], 1);
            subbandOutputY = ifft(blocks.blockYX() .* subbandX + blocks.blockYY() .* subbandY, [], 1);
            clear subbandX subbandY;

            if parityX
                spectrumOutputX = spectrumOutputX + subbandOutputX .* obj.demodulationX_;
                spectrumOutputY = spectrumOutputY + subbandOutputY .* obj.demodulationX_;
            else
                spectrumOutputX = subbandOutputX;
                spectrumOutputY = subbandOutputY;
            end

            clear subbandOutputX subbandOutputY;
        end

        spectrumOutputX = ifft(spectrumOutputX, [], 2);
        spectrumOutputY = ifft(spectrumOutputY, [], 2);

        if parityY
            outputX = outputX + spectrumOutputX .* obj.demodulationY_;
            outputY = outputY + spectrumOutputY .* obj.demodulationY_;
        else
            outputX = spectrumOutputX;
            outputY = spectrumOutputY;
        end

        clear spectrumOutputX spectrumOutputY;
    end

    outputX = outputX(1:rowCountX, 1:columnCountX);
    outputY = outputY(1:rowCountY, 1:columnCountY);

    output = [outputX(:); outputY(:)];
end
