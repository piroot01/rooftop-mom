function spectra = generateBlockwiseSpectra(kernel, generators, faceSizeX, faceSizeY)
% generateBlockwiseSpectra: Returns the block spectra of a circulant approximation of the Green's matrix built block by
%   block: the generators of the 2x2 blocks are formed in space (utils.generateSpatialBlocks) and each is replaced by
%   its circulant approximation. The YX block is the transpose of the XY block, and so is its circulant, which keeps
%   the approximation symmetric for every kernel.
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

    [blockXX, blockXY, ~, blockYY] = rooftopMoM.assembler.preconditioner.utils.generateSpatialBlocks(generators, ...
        faceSizeX, faceSizeY);
    coupling = fft2(kernel(blockXY));

    spectra = rooftopMoM.assembler.greensMatrix.Spectra(fft2(kernel(blockXX)), coupling, transpose_(coupling), ...
        fft2(kernel(blockYY)));
end

function transposed = transpose_(spectrum)
    % The transposed circulant has the generator c(-k, -l), so its spectrum is taken at the negative frequencies.
    if isempty(spectrum)
        transposed = [];
        return;
    end

    transposed = spectrum([1, end:(-1):2], [1, end:(-1):2]);
end
