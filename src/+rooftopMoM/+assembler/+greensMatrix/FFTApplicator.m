classdef FFTApplicator < rooftopMoM.assembler.greensMatrix.Applicator
% FFTApplicator: Applies a block circulant matrix, given by its 2x2 block spectra, by 2D FFTs. Each rooftop component
%   is zero-padded to the spectral grid, transformed, mixed by the blocks, transformed back and cut to its grid.
%   With the spectra of the circulant embedding (utils.generateEmbeddingSpectra) it applies the Green's matrix, with
%   inverted circulant spectra (Spectra.invert) a circulant preconditioner.

    methods
        function obj = FFTApplicator(rooftops, spectra, content)
            obj@rooftopMoM.assembler.greensMatrix.Applicator(rooftops, spectra, content);
        end
    end

    methods
        function output = apply(obj, input)
        % apply: Returns the product with a vector of rooftop coefficients.
        %
        %   Input:
        %       input: Nx1 vector of the rooftop coefficients.
        %
        %   Output:
        %       output: Nx1 vector of the product.

            switch obj.content_
                case rooftopMoM.basisFunctions.RooftopContent.BOTH
                    output = obj.applyBothDirections_(input);
                case rooftopMoM.basisFunctions.RooftopContent.XONLY
                    output = obj.applySingleDirection_(input, 'x', obj.spectra_.blockXX());
                otherwise
                    output = obj.applySingleDirection_(input, 'y', obj.spectra_.blockYY());
            end
        end
    end

    methods (Static)
        function bytes = estimateWorkingMemory(rooftops, sizeX, sizeY, content)
        % estimateWorkingMemory: Returns the bytes of the arrays alive at the peak of one product, the output included,
        %   for the spectra on a sizeX x sizeY grid.

            blockBytes = 16 * sizeX * sizeY;

            if content == rooftopMoM.basisFunctions.RooftopContent.BOTH
                % Both output spectra, the padded spectrum of the y-directed input and the y-directed input itself.
                bytes = 3 * blockBytes + 16 * rooftops.countDir('y');
            else
                % The input spectrum, the product and its inverse transform, partly overlapping.
                bytes = 2.5 * blockBytes;
            end
        end
    end

    methods (Access = private)
        function output = applySingleDirection_(obj, input, direction, block)
        % applySingleDirection_: Applies the matrix of a strip, the rooftops of one direction only.

            [rowCount, columnCount] = obj.rooftops_.gridSize(direction);

            output = ifft2(block .* fft2(reshape(input, rowCount, columnCount), size(block, 1), size(block, 2)));
            output = reshape(output(1:rowCount, 1:columnCount), [], 1);
        end

        function output = applyBothDirections_(obj, input)
        % applyBothDirections_: Applies the matrix of the x- and y-directed rooftops. The arrays are cleared as soon as
        %   possible to keep the peak memory low.

            [rowCountX, columnCountX] = obj.rooftops_.gridSize('x');
            [rowCountY, columnCountY] = obj.rooftops_.gridSize('y');
            [sizeX, sizeY] = obj.spectra_.gridSize();
            countX = rowCountX * columnCountX;

            spectrumX = fft2(reshape(input(1:countX), rowCountX, columnCountX), sizeX, sizeY);
            outputX = obj.spectra_.blockXX() .* spectrumX;
            outputY = obj.spectra_.blockYX() .* spectrumX;
            clear spectrumX;

            spectrumY = fft2(reshape(input((countX + 1):end), rowCountY, columnCountY), sizeX, sizeY);
            outputX = outputX + obj.spectra_.blockXY() .* spectrumY;
            outputY = outputY + obj.spectra_.blockYY() .* spectrumY;
            clear spectrumY;

            outputX = ifft2(outputX);
            outputX = outputX(1:rowCountX, 1:columnCountX);

            outputY = ifft2(outputY);
            outputY = outputY(1:rowCountY, 1:columnCountY);

            output = [outputX(:); outputY(:)];
        end
    end
end
