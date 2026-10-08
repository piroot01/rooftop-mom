function spectralMatrixData = generateSpectralGreensMatrix(magneticXXCirc, magneticYYCirc, electricCirc, rooftopContent, faceSizeX, faceSizeY)
% generateSpectralGreensMatrix: Generates the spectral blocks of the Green's matrix from the circulant embeddings
%   of its parts.
%
%   Inputs:
%       magneticXXCirc: Circulant embedding of the magnetic part of the x-directed rooftops (may be empty).
%       magneticYYCirc: Circulant embedding of the magnetic part of the y-directed rooftops (may be empty).
%       electricCirc: Circulant embedding of the electric part (may be empty).
%       rooftopContent: RooftopContent specifying which blocks are generated.
%       faceSizeX: Face size in the x direction.
%       faceSizeY: Face size in the y direction.
%
%   Output:
%       spectralMatrixData: SpectralMatrixData with the blocks XX, XY, YX and YY.

    if ~isempty(electricCirc)
        [circSizeX, circSizeY] = size(electricCirc);
    elseif ~isempty(magneticXXCirc)
        [circSizeX, circSizeY] = size(magneticXXCirc);
    elseif ~isempty(magneticYYCirc)
        [circSizeX, circSizeY] = size(magneticYYCirc);
    else
        error('At least one of the circulants must be nonempty!');
    end

    [forwardShiftX, ~, forwardShiftY, ~] = ...
        rooftopMoM.core.greensMatrix.utils.generateSpectralShiftVectors(circSizeX, circSizeY, 0, 0);

    spectralMatrixData = rooftopMoM.core.greensMatrix.utils.generateSpectralMatrixData( ...
        fft2(magneticXXCirc), fft2(magneticYYCirc), fft2(electricCirc), ...
        rooftopContent, faceSizeX, faceSizeY, forwardShiftX, forwardShiftY ...
    );
end
