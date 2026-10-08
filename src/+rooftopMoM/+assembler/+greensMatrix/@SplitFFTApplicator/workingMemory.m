function bytes = workingMemory(obj)
% workingMemory: Returns the bytes of the arrays alive at the peak of one product, the output included.

    [fftSizeX, fftSizeY] = obj.spectra_(1).gridSize();
    gridBytes = 16 * fftSizeX * fftSizeY;

    if obj.content_ == rooftopMoM.basisFunctions.RooftopContent.BOTH
        % Both inputs, the first-level spectra and subbands, the mixed x output and the mixing temporary of y.
        bytes = 8 * gridBytes;
    else
        % The input copy, its transform and the output of the even subband.
        bytes = 3 * gridBytes;
    end
end
