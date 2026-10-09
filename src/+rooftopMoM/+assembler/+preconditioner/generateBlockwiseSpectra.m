function spectra = generateBlockwiseSpectra(kernel, generators, faceSizeX, faceSizeY)
% generateBlockwiseSpectra: Returns the block spectra of a circulant approximation of the Green's matrix built block by
%   block: the generators of the four 2x2 blocks are formed in space (utils.generateSpatialBlocks) and each is
%   replaced by its circulant approximation.
%
%   Inputs:
%       kernel: Function handle returning the circulant approximation of a two-sided generator
%           (preconditioner.kernels).
%       generators: rooftopMoM.assembler.greensMatrix.Generators.
%       faceSizeX: Face size in the x direction.
%       faceSizeY: Face size in the y direction.
%
%   Output:
%       spectra: rooftopMoM.assembler.greensMatrix.Spectra on the face grid.

    [blockXX, blockXY, blockYX, blockYY] = rooftopMoM.assembler.preconditioner.utils.generateSpatialBlocks( ...
        generators, faceSizeX, faceSizeY);
    approximate = @(generator) fft2(kernel(rooftopMoM.assembler.preconditioner.utils.extendEvenly(generator)));

    spectra = rooftopMoM.assembler.greensMatrix.Spectra(approximate(blockXX), approximate(blockXY), ...
        approximate(blockYX), approximate(blockYY));
end
