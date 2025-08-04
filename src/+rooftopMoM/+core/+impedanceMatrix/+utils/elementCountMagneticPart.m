function count = elementCountMagneticPart(context, remainingMemory)
% elementCountMagneticPart: Computes the element count for the magnetic part of the impedance matrix, for use by the assembler chunk planner.
%
%   Inputs:
%       context: AssemblerContext instance.
%       remainingMemory: Memory instance specifying the memory available to the internal chunk-loop of the magnetic assembler.
%
%   Output:
%       count: Number of elements to assign to the chunk planner.

    % The most efficient is to pick square chunk size, hence we assume that a = b.
    % The constant in front of the ab = a^2 consists of:
    % rowRooftopPoints: D * (a * q^2) * 3 * 2
    % columnRooftopPoints: D * (b * q^2) * 3 * 2
    % rowQuadraturePoints: D * (a * q^2) * 3 * 2
    % columnQuadraturePoints: D * (b * q^2) * 3 * 2
    % distanceMatrix: D * (a * b * q^4) * 4
    % integralKernel: CD * (a * b * q^4) * 4
    % idMask: UI32 * a * b
    % sizeMask: UI32 * a * b
    % singularityMask: L * a * b

    order = context.quadrature().order();

    complexMemory = rooftopMoM.DataTypeMemory.memory(rooftopMoM.DataTypeMemory.ComplexDouble);
    doubleMemory = rooftopMoM.DataTypeMemory.memory(rooftopMoM.DataTypeMemory.Double);
    uint32Memory = rooftopMoM.DataTypeMemory.memory(rooftopMoM.DataTypeMemory.UInt32);
    logicalMemory = rooftopMoM.DataTypeMemory.memory(rooftopMoM.DataTypeMemory.Logical);

    % All terms with a linear form of a/b
    c1 = doubleMemory * (order^2 * 3 * 2 + order^2 * 3 * 2 + order^2 * 3 * 2 + order^2 * 3 * 2);

    % All terms with a quadratic form a * b
    c2 = doubleMemory * (order^4 * 4) + complexMemory * (order^4 * 4) + uint32Memory * 2 + logicalMemory;

    c1Bytes = c1.bytes();
    c2Bytes = c2.bytes();

    % Next we solve c1 * a + c2 * a^2 = remainingMemory
    denominator = sqrt(c1Bytes^2 + 4 * c2Bytes * remainingMemory.bytes()) - c1Bytes;
    numerator = 2 * c2Bytes;

    count = floor((denominator / numerator)^2);
end
