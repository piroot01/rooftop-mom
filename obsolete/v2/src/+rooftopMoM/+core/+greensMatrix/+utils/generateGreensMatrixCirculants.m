function [magneticXXCirc, magneticYYCirc, electricCirc] = generateGreensMatrixCirculants(matrixData, cX, cY)
    if nargin < 2
        [cX, cY] = size(matrixData.electric());
    end

    switch matrixData.rooftopContent()
        case rooftopMoM.basisFunctions.RooftopContent.XONLY
            magneticXXCirc = rooftopMoM.core.greensMatrix.utils.embedToeplitz4ToCirculant(matrixData.magneticXX(), cX, cY);
            magneticYYCirc = [];
            electricCirc = rooftopMoM.core.greensMatrix.utils.embedToeplitz4ToCirculant(matrixData.electric(), cX, cY);
        case rooftopMoM.basisFunctions.RooftopContent.YONLY
            magneticXXCirc = [];
            magneticYYCirc = rooftopMoM.core.greensMatrix.utils.embedToeplitz4ToCirculant(matrixData.magneticYY(), cX, cY);
            electricCirc = rooftopMoM.core.greensMatrix.utils.embedToeplitz4ToCirculant(matrixData.electric(), cX, cY);
        case rooftopMoM.basisFunctions.RooftopContent.BOTH
            magneticXXCirc = rooftopMoM.core.greensMatrix.utils.embedToeplitz4ToCirculant(matrixData.magneticXX(), cX, cY);
            magneticYYCirc = rooftopMoM.core.greensMatrix.utils.embedToeplitz4ToCirculant(matrixData.magneticYY(), cX, cY);
            electricCirc = rooftopMoM.core.greensMatrix.utils.embedToeplitz4ToCirculant(matrixData.electric(), cX, cY);
        otherwise
            error('Unknown rooftop content: %s', rooftopContent);
    end
end
