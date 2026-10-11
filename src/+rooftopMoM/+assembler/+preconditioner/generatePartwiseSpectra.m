function spectra = generatePartwiseSpectra(kernel, generators, faceSizeX, faceSizeY)
% generatePartwiseSpectra: Returns the block spectra of a circulant approximation of the Green's matrix built part by
%   part: each generator, zero-padded to the face grid, is replaced by its circulant approximation, and the blocks
%   follow as for the Green's matrix itself.
%
%   Inputs:
%       kernel: Function handle returning the circulant approximation of a two-sided generator
%           (preconditioner.kernels); a kernel with a second input gets the one-sided generators and true.
%       generators: rooftopMoM.assembler.greensMatrix.Generators.
%       faceSizeX: Face size in the x direction.
%       faceSizeY: Face size in the y direction.
%
%   Output:
%       spectra: rooftopMoM.assembler.greensMatrix.Spectra on the face grid.

    [countX, countY] = size(generators.electric());
    pad = @(generator) rooftopMoM.assembler.utils.padWithZeros(generator, countX, countY);
    approximate = @(generator) fft2(kernel(rooftopMoM.assembler.preconditioner.utils.extendEvenly(pad(generator))));

    % The generators are even, so the kernels that take a one-sided even generator as well, T. Chan and R. Chan, fold
    % them without their two-sided extension.
    if nargin(kernel) > 1
        approximate = @(generator) fft2(kernel(pad(generator), true));
    end

    [shiftX, shiftY] = rooftopMoM.assembler.greensMatrix.utils.generateShiftPhases(countX, countY, 0, 0);

    spectra = rooftopMoM.assembler.greensMatrix.utils.generateBlockSpectra(approximate(generators.magneticXX()), ...
        approximate(generators.magneticYY()), approximate(generators.electric()), generators.content(), ...
        faceSizeX, faceSizeY, shiftX, shiftY);
end
