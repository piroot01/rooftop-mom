function generator = generateElectricPart(context)
% generateElectricPart: Returns the Toeplitz generator of the scalar potential part of the Green's matrix: the
%   interactions of every face with the face at the origin.
%
%   Input:
%       context: rooftopMoM.assembler.AssemblerContext.
%
%   Output:
%       generator: countX x countY generator.

    primitive = context.rooftops().primitive();
    weights = context.rooftops().quadrature().weights4D();
    staticIntegral = rooftopMoM.assembler.singularIntegrals.type3(primitive.sizeX(), primitive.sizeY());

    generator = rooftopMoM.assembler.greensMatrix.utils.integrateFacePairs(context, ...
        [primitive.countX(), primitive.countY()], [0; 0], weights, staticIntegral);

    faceArea = primitive.sizeX() * primitive.sizeY();
    generator = 1i * rooftopMoM.assembler.Constants.Z0 / (4 * pi * context.wavenumber()) * faceArea^2 * generator;
end
