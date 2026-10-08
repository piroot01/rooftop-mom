classdef LeanSplitFFTApplicator < rooftopMoM.assembler.greensMatrix.SplitFFTApplicator
% LeanSplitFFTApplicator: Applies the Green's matrix by the split FFT with the least memory.
%   In 2D every contribution (subband, input component, output component) is transformed alone in one grid-sized
%   buffer and added to the output. The peak is the output and a single grid array, a quarter of the working memory
%   of the FFTApplicator, at the price of 64 FFTs instead of the 24 of the SplitFFTApplicator. Strips use the
%   inherited product, which is already memory minimal.

    properties (Constant, Access = private)
        % Number of column chunks copied between the grid buffer and the vectors; it bounds the temporaries.
        chunkCount_ = 32;
    end

    methods
        function obj = LeanSplitFFTApplicator(rooftops, generators)
            obj@rooftopMoM.assembler.greensMatrix.SplitFFTApplicator(rooftops, generators);
        end
    end

    methods
        bytes = workingMemory(obj);
    end

    methods (Access = protected)
        output = applyBothDirections_(obj, input);
    end

    methods (Access = private)
        [rowCounts, columnCounts, offsets] = componentLayout_(obj);
    end
end
