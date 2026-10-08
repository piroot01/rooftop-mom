function spectralMatrixData = generateSpectralMatrixData(magneticXXSpectr, magneticYYSpectr, electricSpectr, rooftopContent, faceSizeX, faceSizeY, forwardShiftX, forwardShiftY)
% generateSpectralMatrixData: Assembles the spectral blocks of the Green's matrix from the spectra of its parts.
%
%   Inputs:
%       magneticXXSpectr: Spectrum of the magnetic part of the x-directed rooftops (empty is treated as zero).
%       magneticYYSpectr: Spectrum of the magnetic part of the y-directed rooftops (empty is treated as zero).
%       electricSpectr: Spectrum of the electric part (empty is treated as zero).
%       rooftopContent: RooftopContent specifying which blocks are assembled.
%       faceSizeX: Face size in the x direction.
%       faceSizeY: Face size in the y direction.
%       forwardShiftX: Column vector of the forward shift in the x direction sampled at the spectrum frequencies.
%       forwardShiftY: Row vector of the forward shift in the y direction sampled at the spectrum frequencies.
%
%   Output:
%       spectralMatrixData: SpectralMatrixData with the blocks XX, XY, YX and YY (the inactive ones are empty).

    if isempty(magneticXXSpectr)
        magneticXXSpectr = 0;
    end

    if isempty(magneticYYSpectr)
        magneticYYSpectr = 0;
    end

    if isempty(electricSpectr)
        electricSpectr = 0;
    end

    isXActive = (rooftopContent == rooftopMoM.basisFunctions.RooftopContent.XONLY) || ...
        (rooftopContent == rooftopMoM.basisFunctions.RooftopContent.BOTH);

    isYActive = (rooftopContent == rooftopMoM.basisFunctions.RooftopContent.YONLY) || ...
        (rooftopContent == rooftopMoM.basisFunctions.RooftopContent.BOTH);

    blockXX = []; blockXY = []; blockYX = []; blockYY = [];

    % NOTE: (1 - f) .* (conj(f) - 1) = -|1 - f|^2 keeps the diagonal blocks real for real parts.
    if isXActive
        blockXX = magneticXXSpectr - (abs(1 - forwardShiftX) / faceSizeX).^2 .* electricSpectr;
    end

    if isYActive
        blockYY = magneticYYSpectr - (abs(1 - forwardShiftY) / faceSizeY).^2 .* electricSpectr;
    end

    if isXActive && isYActive
        blockXY = ((1 - forwardShiftX) / (faceSizeX * faceSizeY)) .* (conj(forwardShiftY) - 1) .* electricSpectr;
        blockYX = ((conj(forwardShiftX) - 1) / (faceSizeX * faceSizeY)) .* (1 - forwardShiftY) .* electricSpectr;
    end

    spectralMatrixData = rooftopMoM.core.greensMatrix.SpectralMatrixData( ...
        'blockXX', blockXX, ...
        'blockXY', blockXY, ...
        'blockYX', blockYX, ...
        'blockYY', blockYY ...
    );
end
