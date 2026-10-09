function spectra = generateHybridSpectra(kernel, generators, faceSizeX, faceSizeY)
% generateHybridSpectra: Returns the block spectra of a circulant approximation of the Green's matrix that combines
%   the two builders: the diagonal blocks are approximated in space as in generateBlockwiseSpectra, the coupling
%   blocks follow from the circulant approximation of the electric generator as in generatePartwiseSpectra.
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

    [blockXX, ~, ~, blockYY] = rooftopMoM.assembler.preconditioner.utils.generateSpatialBlocks(generators, ...
        faceSizeX, faceSizeY);
    approximate = @(generator) fft2(kernel(rooftopMoM.assembler.preconditioner.utils.extendEvenly(generator)));

    [countX, countY] = size(generators.electric());
    [shiftX, shiftY] = rooftopMoM.assembler.greensMatrix.utils.generateShiftPhases(countX, countY, 0, 0);
    coupling = rooftopMoM.assembler.greensMatrix.utils.generateBlockSpectra([], [], ...
        approximate(generators.electric()), generators.content(), faceSizeX, faceSizeY, shiftX, shiftY);

    spectra = rooftopMoM.assembler.greensMatrix.Spectra(approximate(blockXX), coupling.blockXY(), ...
        coupling.blockYX(), approximate(blockYY));
end
