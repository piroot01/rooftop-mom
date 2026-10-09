classdef LeanSplitFFTApplicator < rooftopMoM.assembler.greensMatrix.SplitFFTApplicator
% LeanSplitFFTApplicator: Applies the Green's matrix by the split FFT with the least memory, from the subband spectra
%   of utils.generateSubbandSpectra.
%   In 2D every contribution (subband, input component, output component) is transformed alone in one grid-sized
%   buffer and added to the output. The peak is the output and a single grid array, a quarter of the working memory
%   of the FFTApplicator, at the price of 64 FFTs instead of the 24 of the SplitFFTApplicator. Strips use the
%   inherited product, which is already memory minimal.

    properties (Constant, Access = private)
        % Number of column chunks copied between the grid buffer and the vectors; it bounds the temporaries.
        chunkCount_ = 32;
    end

    methods
        function obj = LeanSplitFFTApplicator(rooftops, spectra, content)
            obj@rooftopMoM.assembler.greensMatrix.SplitFFTApplicator(rooftops, spectra, content);
        end
    end

    methods (Static)
        bytes = estimateWorkingMemory(rooftops, sizeX, sizeY, content);
    end

    methods (Access = protected)
        output = applyBothDirections_(obj, input);
    end

    methods (Static, Access = private)
        [rowCounts, columnCounts, offsets] = componentLayout_(rooftops);
    end
end
