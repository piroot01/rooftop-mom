function bytes = workingMemory(obj)
% workingMemory: Returns the bytes of the arrays alive at the peak of one product, the output included.

    if obj.content_ ~= rooftopMoM.basisFunctions.RooftopContent.BOTH
        bytes = workingMemory@rooftopMoM.assembler.greensMatrix.SplitFFTApplicator(obj);
        return;
    end

    [fftSizeX, fftSizeY] = obj.spectra_(1).gridSize();
    [rowCounts, columnCounts] = obj.componentLayout_();

    % The largest column chunk is copied up to three times at once (slice, output range and their sum).
    chunkElements = max(rowCounts .* ceil(columnCounts ./ min(obj.chunkCount_, columnCounts)));

    bytes = 16 * (obj.rooftops_.count() + fftSizeX * fftSizeY + 3 * chunkElements);
end
