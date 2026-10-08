function generator = generateMagneticPart(context, direction)
% generateMagneticPart: Returns the Toeplitz generator of the vector potential part of the Green's matrix for the
%   rooftops of one direction: the interactions of every rooftop with the reference rooftop at the origin.
%
%   Inputs:
%       context: rooftopMoM.assembler.AssemblerContext.
%       direction: 'x' or 'y'.
%
%   Output:
%       generator: Generator over the rooftop grid of the direction (UniformRooftops.gridSize).

    rooftops = context.rooftops();
    primitive = rooftops.primitive();
    rooftop = rooftops.referenceRooftop(direction);
    weights = rooftops.quadrature().weights2D();
    [rowCount, columnCount] = rooftops.gridSize(direction);

    % Face 1 of a rooftop is the rising half, face 2 the falling one. On a shared face, equal halves give the static
    % integral type1 and different halves type2.
    faceShifts = [zeros(2, 1), rooftop.secondFaceOffset];
    sameHalves = rooftopMoM.assembler.singularIntegrals.type1(rooftop.length, rooftop.height);
    differentHalves = rooftopMoM.assembler.singularIntegrals.type2(rooftop.length, rooftop.height);
    staticIntegrals = [sameHalves, differentHalves; differentHalves, sameHalves];

    generator = complex(zeros(rowCount, columnCount));

    for face = 1:2
        for referenceFace = 1:2
            pairWeights = (weights .* rooftop.amplitudes(:, face)) * (weights .* rooftop.amplitudes(:, referenceFace)).';

            generator = generator + rooftopMoM.assembler.greensMatrix.utils.integrateFacePairs(context, ...
                [rowCount, columnCount], faceShifts(:, face) - faceShifts(:, referenceFace), pairWeights, ...
                staticIntegrals(face, referenceFace));
        end
    end

    faceArea = primitive.sizeX() * primitive.sizeY();
    generator = 1i * context.wavenumber() * rooftopMoM.assembler.Constants.Z0 / (4 * pi) * faceArea^2 * generator;
end
