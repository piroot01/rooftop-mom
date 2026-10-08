function spectra = generateEmbeddingSpectra(generators, faceSizeX, faceSizeY)
% generateEmbeddingSpectra: Returns the block spectra of the Green's matrix on the grid of its circulant embedding,
%   for the FFTApplicator. In every direction with more than one face, the symmetric Toeplitz generator
%   t(0), ..., t(n - 1) is mirrored into the circulant t(0), ..., t(n - 1), t(n - 2), ..., t(1) of the period 2n - 2.
%   The magnetic generators are zero-padded to the size of the electric one first, so all blocks share one grid.
%
%   Inputs:
%       generators: rooftopMoM.assembler.greensMatrix.Generators.
%       faceSizeX: Face size in the x direction.
%       faceSizeY: Face size in the y direction.
%
%   Output:
%       spectra: rooftopMoM.assembler.greensMatrix.Spectra.

    start = tic;
    [countX, countY] = size(generators.electric());
    mirroredX = [1:countX, (countX - 1):(-1):2];
    mirroredY = [1:countY, (countY - 1):(-1):2];
    embed = @(generator) fft2(mirror_(rooftopMoM.assembler.utils.padWithZeros(generator, countX, countY), ...
        mirroredX, mirroredY));

    [shiftX, shiftY] = rooftopMoM.assembler.greensMatrix.utils.generateShiftPhases( ...
        numel(mirroredX), numel(mirroredY), 0, 0);

    spectra = rooftopMoM.assembler.greensMatrix.utils.generateBlockSpectra(embed(generators.magneticXX()), ...
        embed(generators.magneticYY()), embed(generators.electric()), generators.content(), faceSizeX, faceSizeY, ...
        shiftX, shiftY);

    rooftopMoM.logging.Logger.instance().info('Generated the embedding spectra in %.3f s (%.1f MB).', toc(start), ...
        spectra.memory() / 2^20);
end

function mirrored = mirror_(generator, mirroredX, mirroredY)
    if isempty(generator)
        mirrored = [];
        return;
    end

    mirrored = generator(mirroredX, mirroredY);
end
