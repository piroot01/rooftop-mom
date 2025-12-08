function determinant = generateSpectralGreensMatrixDeterminant(gxxf, gxyf, gyxf, gyyf, rooftopContent)
    switch rooftopContent
        case rooftopMoM.basisFunctions.RooftopContent.XONLY
            determinant = [];
        case rooftopMoM.basisFunctions.RooftopContent.YONLY
            determinant = [];
        case rooftopMoM.basisFunctions.RooftopContent.BOTH
            determinant = gxxf .* gyyf - gxyf .* gyxf;
        otherwise
            error('Unknown rooftop content: %s', rooftopContent);
    end
end
