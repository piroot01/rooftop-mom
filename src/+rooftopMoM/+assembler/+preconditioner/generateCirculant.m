function preconditioner = generateCirculant(rooftops, generators, kernel, builder)
% generateCirculant: Returns the circulant preconditioner of the Green's matrix: the FFT product with the inverse of a
%   block circulant approximation.
%
%   Inputs:
%       rooftops: rooftopMoM.basisFunctions.UniformRooftops.
%       generators: rooftopMoM.assembler.greensMatrix.Generators.
%       kernel: Function handle of the circulant approximation of a generator, e.g.
%           @rooftopMoM.assembler.preconditioner.kernels.generateStrang.
%       builder: Function handle of the spectra builder, e.g.
%           @rooftopMoM.assembler.preconditioner.generatePartwiseSpectra.
%
%   Output:
%       preconditioner: rooftopMoM.assembler.greensMatrix.FFTApplicator of the inverted spectra.

    start = tic;
    primitive = rooftops.primitive();
    spectra = builder(kernel, generators, primitive.sizeX(), primitive.sizeY());

    preconditioner = rooftopMoM.assembler.greensMatrix.FFTApplicator(rooftops, spectra.invert(generators.content()), ...
        generators.content());

    rooftopMoM.logging.Logger.instance().info('Built the circulant preconditioner in %.3f s (%.1f MB).', ...
        toc(start), preconditioner.memory() / 2^20);
end
