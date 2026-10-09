function [stored, working] = estimatePreconditioner(rooftops)
% estimatePreconditioner: Returns the memory of a circulant preconditioner before it is built: the bytes of its
%   spectra on the face grid and the bytes alive at the peak of one application.
%
%   Input:
%       rooftops: rooftopMoM.basisFunctions.UniformRooftops.
%
%   Outputs:
%       stored: Bytes of the spectra.
%       working: Bytes of the arrays alive at the peak of one application.

    primitive = rooftops.primitive();
    content = rooftops.content();
    blockCount = 1 + 3 * (content == rooftopMoM.basisFunctions.RooftopContent.BOTH);

    stored = 16 * blockCount * primitive.count();
    working = rooftopMoM.assembler.greensMatrix.FFTApplicator.estimateWorkingMemory(rooftops, primitive.countX(), ...
        primitive.countY(), content);
end
