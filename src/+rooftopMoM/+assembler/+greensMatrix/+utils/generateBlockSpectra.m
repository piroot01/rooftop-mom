function spectra = generateBlockSpectra(magneticXX, magneticYY, electric, content, faceSizeX, faceSizeY, shiftX, shiftY)
% generateBlockSpectra: Combines the spectra of the three generators into the 2x2 block spectra of the Green's
%   matrix. The scalar potential acts through the divergence of the rooftops, which is the difference
%   (1 - shift) / faceSize in the spectral domain.
%
%   Inputs:
%       magneticXX: Spectrum of the magnetic part of the x-directed rooftops (empty is treated as zero).
%       magneticYY: Spectrum of the magnetic part of the y-directed rooftops (empty is treated as zero).
%       electric: Spectrum of the electric part.
%       content: rooftopMoM.basisFunctions.RooftopContent, which blocks are formed.
%       faceSizeX: Face size in the x direction.
%       faceSizeY: Face size in the y direction.
%       shiftX: Column vector of the phase of a unit shift in x at the frequencies of the grid.
%       shiftY: Row vector of the phase of a unit shift in y at the frequencies of the grid.
%
%   Output:
%       spectra: rooftopMoM.assembler.greensMatrix.Spectra, the blocks of a missing direction empty.

    if isempty(magneticXX)
        magneticXX = 0;
    end

    if isempty(magneticYY)
        magneticYY = 0;
    end

    hasX = content ~= rooftopMoM.basisFunctions.RooftopContent.YONLY;
    hasY = content ~= rooftopMoM.basisFunctions.RooftopContent.XONLY;

    blockXX = [];
    blockXY = [];
    blockYX = [];
    blockYY = [];

    % (1 - shift) .* (conj(shift) - 1) = -|1 - shift|^2 keeps the diagonal blocks real for real generators.
    if hasX
        blockXX = magneticXX - (abs(1 - shiftX) / faceSizeX).^2 .* electric;
    end

    if hasY
        blockYY = magneticYY - (abs(1 - shiftY) / faceSizeY).^2 .* electric;
    end

    if hasX && hasY
        blockXY = ((1 - shiftX) / (faceSizeX * faceSizeY)) .* (conj(shiftY) - 1) .* electric;
        blockYX = ((conj(shiftX) - 1) / (faceSizeX * faceSizeY)) .* (1 - shiftY) .* electric;
    end

    spectra = rooftopMoM.assembler.greensMatrix.Spectra(blockXX, blockXY, blockYX, blockYY);
end
