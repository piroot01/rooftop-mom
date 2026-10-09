function bytes = estimateWorkingMemory(~, sizeX, sizeY, content)
% estimateWorkingMemory: Returns the bytes of the arrays alive at the peak of one product, the output included, for
%   the subband spectra on a sizeX x sizeY grid.
%
%   Inputs:
%       ~: The rooftops, unused; the peak depends on the grid only.
%       sizeX, sizeY: Size of the fast FFT grid of the subbands.
%       content: rooftopMoM.basisFunctions.RooftopContent.
%
%   Output:
%       bytes: Peak memory of one product.

    gridBytes = 16 * sizeX * sizeY;

    if content == rooftopMoM.basisFunctions.RooftopContent.BOTH
        % Both inputs, the first-level spectra and subbands, the mixed x output and the mixing temporary of y.
        bytes = 8 * gridBytes;
    else
        % The input copy, its transform and the output of the even subband.
        bytes = 3 * gridBytes;
    end
end
