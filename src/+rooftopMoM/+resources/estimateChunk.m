function bytes = estimateChunk(rooftops, orders)
% estimateChunk: Returns the bytes of the largest chunk the assembly of the Green's matrix can use, before it is
%   assembled: its near or its far pass in one chunk. A larger chunk changes nothing.
%
%   Inputs:
%       rooftops: rooftopMoM.basisFunctions.UniformRooftops.
%       orders: 1x2 Gauss-Legendre orders [near, far] of the cell pairs.
%
%   Output:
%       bytes: Memory of the chunk.

    bytes = rooftopMoM.assembler.greensMatrix.MomentAssembler.estimateWorkingMemory(rooftops, orders);
end
