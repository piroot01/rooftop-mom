function bytes = estimateWorkingMemory(rooftops, sizeX, sizeY, content)
% estimateWorkingMemory: Returns the bytes of the arrays alive at the peak of one product, the output included, for
%   the subband spectra on a sizeX x sizeY grid.
%
%   Inputs:
%       rooftops: rooftopMoM.basisFunctions.UniformRooftops.
%       sizeX, sizeY: Size of the fast FFT grid of the subbands.
%       content: rooftopMoM.basisFunctions.RooftopContent.
%
%   Output:
%       bytes: Peak memory of one product.

    if content ~= rooftopMoM.basisFunctions.RooftopContent.BOTH
        bytes = rooftopMoM.assembler.greensMatrix.SplitFFTApplicator.estimateWorkingMemory(rooftops, sizeX, sizeY, ...
            content);
        return;
    end

    [rowCounts, columnCounts] = rooftopMoM.assembler.greensMatrix.LeanSplitFFTApplicator.componentLayout_(rooftops);
    chunkCount = rooftopMoM.assembler.greensMatrix.LeanSplitFFTApplicator.chunkCount_;

    % The largest column chunk is copied up to three times at once (slice, output range and their sum).
    chunkElements = max(rowCounts .* ceil(columnCounts ./ min(chunkCount, columnCounts)));

    bytes = 16 * (rooftops.count() + sizeX * sizeY + 3 * chunkElements);
end
