classdef SplitFFTApplicator < rooftopMoM.assembler.greensMatrix.Applicator
% SplitFFTApplicator: Applies the Green's matrix by the split (subband) FFT with the fewest operations, from the
%   subband spectra of utils.generateSubbandSpectra.
%   Only arrays of the (fast FFT) grid size are transformed, with no zero-padding to the doubled circulant. In 2D the
%   transforms of the first split direction are shared by two subbands (24 FFTs in total), which keeps more arrays
%   alive than the LeanSplitFFTApplicator.

    properties (Access = protected)
        % Odd subband modulations in the x (column) and the y (row) direction, and their complex conjugates.
        modulationX_;
        modulationY_;
        demodulationX_;
        demodulationY_;
    end

    methods
        function obj = SplitFFTApplicator(rooftops, spectra, content)
            obj@rooftopMoM.assembler.greensMatrix.Applicator(rooftops, spectra, content);

            [fftSizeX, fftSizeY] = spectra(1).gridSize();
            obj.modulationX_ = rooftopMoM.assembler.greensMatrix.utils.generateSubbandModulation(fftSizeX, 1);
            obj.modulationY_ = rooftopMoM.assembler.greensMatrix.utils.generateSubbandModulation(fftSizeY, 2);
            obj.demodulationX_ = conj(obj.modulationX_);
            obj.demodulationY_ = conj(obj.modulationY_);
        end
    end

    methods
        output = apply(obj, input);
    end

    methods (Static)
        bytes = estimateWorkingMemory(rooftops, sizeX, sizeY, content);
    end

    methods (Access = protected)
        output = applySingleDirection_(obj, input, evenBlock, oddBlock, modulation, demodulation);

        output = applyBothDirections_(obj, input);
    end
end
