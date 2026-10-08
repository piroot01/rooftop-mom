function spectra = generateBlockwiseSpectra(kernel, generators, faceSizeX, faceSizeY)
% generateBlockwiseSpectra: Returns the block spectra of a circulant approximation of the Green's matrix built block by
%   block: the generators of the four 2x2 blocks are formed in space (utils.generateSpatialBlocks) and each is
%   replaced by its circulant approximation.
%
%   Inputs:
%       kernel: Function handle returning the circulant approximation of a generator (preconditioner.kernels).
%       generators: rooftopMoM.assembler.greensMatrix.Generators.
%       faceSizeX: Face size in the x direction.
%       faceSizeY: Face size in the y direction.
%
%   Output:
%       spectra: rooftopMoM.assembler.greensMatrix.Spectra on the face grid.

    [blockXX, blockXY, blockYX, blockYY] = rooftopMoM.assembler.preconditioner.utils.generateSpatialBlocks( ...
        generators, faceSizeX, faceSizeY);

    spectra = rooftopMoM.assembler.greensMatrix.Spectra(fft2(kernel(blockXX)), fft2(kernel(blockXY)), ...
        fft2(kernel(blockYX)), fft2(kernel(blockYY)));
end
