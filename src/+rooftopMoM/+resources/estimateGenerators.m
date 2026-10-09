function bytes = estimateGenerators(rooftops)
% estimateGenerators: Returns the bytes of the generators of the Green's matrix before they are assembled, one complex
%   value per rooftop and per face.
%
%   Input:
%       rooftops: rooftopMoM.basisFunctions.UniformRooftops.
%
%   Output:
%       bytes: Memory of the generators.

    bytes = 16 * (rooftops.count() + rooftops.primitive().count());
end
