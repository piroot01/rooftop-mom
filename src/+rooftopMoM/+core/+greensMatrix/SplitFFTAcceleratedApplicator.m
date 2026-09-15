classdef SplitFFTAcceleratedApplicator < rooftopMoM.core.greensMatrix.GreensMatrixApplicator
% SplitFFTAcceleratedApplicator: Applies the Green's matrix using the split (subband) FFT with the fewest operations.
%   It works with the subband spectra from MatrixData.generateSplitSpectralGreensMatrix, so only arrays of the (fast
%   FFT) grid size are transformed and no zero-padding to the doubled circulant is needed. In 2D the transforms of
%   the first split direction are shared by two subbands (24 FFTs in total), which keeps more arrays alive than the
%   MemoryEfficientSplitFFTAcceleratedApplicator.

    properties (Access = protected)
        % Odd subband modulations in the x (column) and the y (row) direction and their complex conjugates.
        modulationX_;
        modulationY_;
        demodulationX_;
        demodulationY_;
    end

    methods
        function obj = SplitFFTAcceleratedApplicator(context, spectralMatrixData, rooftopContent)
            obj@rooftopMoM.core.greensMatrix.GreensMatrixApplicator(context, spectralMatrixData, rooftopContent);
        end
    end

    methods
        function generate(obj)
            [fftSizeX, fftSizeY] = obj.spectralMatrixData_(1).blockSize();

            obj.modulationX_ = rooftopMoM.core.greensMatrix.utils.generateSubbandModulation(fftSizeX, 1);
            obj.modulationY_ = rooftopMoM.core.greensMatrix.utils.generateSubbandModulation(fftSizeY, 2);

            obj.demodulationX_ = conj(obj.modulationX_);
            obj.demodulationY_ = conj(obj.modulationY_);
        end

        function output = apply(obj, input)
            switch obj.rooftopContent_
                case rooftopMoM.basisFunctions.RooftopContent.XONLY
                    output = obj.applySingleDir_(input, ...
                        obj.spectralMatrixData_(1).blockXX(), obj.spectralMatrixData_(2).blockXX(), ...
                        obj.modulationX_, obj.demodulationX_);
                case rooftopMoM.basisFunctions.RooftopContent.YONLY
                    output = obj.applySingleDir_(input, ...
                        obj.spectralMatrixData_(1).blockYY(), obj.spectralMatrixData_(2).blockYY(), ...
                        obj.modulationY_, obj.demodulationY_);
                case rooftopMoM.basisFunctions.RooftopContent.BOTH
                    output = obj.applyBothDirs_(input);
                otherwise
                    error('Unknown rooftop content: %s', char(obj.rooftopContent_));
            end
        end

        function memory = workingMemory(obj)
        % workingMemory: Returns the peak memory of the arrays alive during one application (output included).
        %
        %   Output:
        %       memory: rooftopMoM.Memory of the peak working arrays.

            [fftSizeX, fftSizeY] = obj.spectralMatrixData_(1).blockSize();
            gridBytes = fftSizeX * fftSizeY * ...
                rooftopMoM.DataTypeMemory.memory(rooftopMoM.DataTypeMemory.ComplexDouble).bytes();

            if obj.rooftopContent_ == rooftopMoM.basisFunctions.RooftopContent.BOTH
                % Both inputs, level-1 spectra and subbands, the mixed x output and the mixing temporary of y.
                memory = rooftopMoM.Memory(8 * gridBytes, 'B');
            else
                % Input copy, its transform and the output of the even subband.
                memory = rooftopMoM.Memory(3 * gridBytes, 'B');
            end
        end
    end

    methods (Access = protected)
        function output = applySingleDir_(~, input, evenBlock, oddBlock, modulation, demodulation)
        % applySingleDir_: Applies the Green's matrix of the rooftops oriented in a single direction.
        %
        %   Inputs:
        %       input: Nx1 vector of the rooftop coefficients.
        %       evenBlock: Spectral block of the even subband (column vector for XONLY, row vector for YONLY).
        %       oddBlock: Spectral block of the odd subband with the same shape as evenBlock.
        %       modulation: Odd subband modulation with the same shape as evenBlock.
        %       demodulation: Complex conjugate of the modulation.
        %
        %   Output:
        %       output: Nx1 vector of the Green's matrix product.

            count = numel(input);

            % Zero-extend the input to the subband size and orient it as the blocks.
            input(numel(evenBlock)) = 0;
            input = reshape(input, size(evenBlock));

            output = ifft(evenBlock .* fft(input));

            % The odd subband is transformed in the input buffer, which is not needed afterwards.
            input = ifft(oddBlock .* fft(input .* modulation)) .* demodulation;
            output = output + input;
            clear input;

            output = reshape(output(1:count), [], 1);
        end

        function output = applyBothDirs_(obj, input)
        % applyBothDirs_: Applies the Green's matrix of the x- and y-directed rooftops.
        %
        %   Input:
        %       input: Nx1 vector of the rooftop coefficients, the x-directed rooftops first.
        %
        %   Output:
        %       output: Nx1 vector of the Green's matrix product.
        %
        %   NOTE: The y direction is split first (2 FFTs per component) and the x direction second (4 FFTs per
        %       component along the contiguous columns), which is the faster order. Both components are written
        %       out explicitly and the arrays are cleared as soon as possible to keep the operations in place.

            primitive = obj.context_.rooftops().primitive();
            countX = double(primitive.countX());
            countY = double(primitive.countY());

            [fftSizeX, fftSizeY] = obj.spectralMatrixData_(1).blockSize();
            inputCountX = (countX - 1) * countY;

            % Zero-extend both components to the subband grid.
            inputX = reshape(input(1:inputCountX), countX - 1, countY);
            inputX(fftSizeX, fftSizeY) = 0;

            inputY = reshape(input((inputCountX + 1):end), countX, countY - 1);
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

                    blocks = obj.spectralMatrixData_(1 + parityX + 2 * parityY);

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

            outputX = outputX(1:(countX - 1), 1:countY);
            outputY = outputY(1:countX, 1:(countY - 1));

            output = [outputX(:); outputY(:)];
        end
    end
end
