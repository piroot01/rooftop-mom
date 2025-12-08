function [gxxf, gxyf, gyxf, gyyf] = generateSpectralGreensMatrix(magneticXX, magneticYY, electric, rooftopContent, faceSizeX, faceSizeY)
    [sizeX, sizeY] = size(electric);

    switch rooftopContent
        case rooftopMoM.basisFunctions.RooftopContent.XONLY
            magneticXXSpectr = fft2(magneticXX);
            electricSpectr = fft2(electric);

            [forwardShiftX, backwardShiftX, ~, ~] = ...
                rooftopMoM.core.greensMatrix.utils.generateSpectralShiftOperators(sizeX, sizeY);

            gxxf = magneticXXSpectr + (1 / faceSizeX)^2 .* (1 - forwardShiftX) .* (backwardShiftX - 1) .* electricSpectr;
            gxyf = []; gyxf = []; gyyf = [];
        case rooftopMoM.basisFunctions.RooftopContent.YONLY
            magneticYYSpectr = fft2(magneticYY);
            electricSpectr = fft2(electric);

            [~, ~, forwardShiftY, backwardShiftY] = ...
                rooftopMoM.core.greensMatrix.utils.generateSpectralShiftOperators(sizeX, sizeY);

            gyyf = magneticYYSpectr + (1 / faceSizeY)^2 .* (1 - forwardShiftY) .* (backwardShiftY - 1) .* electricSpectr;
            gxxf = []; gxyf = []; gyxf = [];
        case rooftopMoM.basisFunctions.RooftopContent.BOTH
            magneticXXSpectr = fft2(magneticXX);
            magneticYYSpectr = fft2(magneticYY);
            electricSpectr = fft2(electric);

            [forwardShiftX, backwardShiftX, forwardShiftY, backwardShiftY] = ...
                rooftopMoM.core.greensMatrix.utils.generateSpectralShiftOperators(sizeX, sizeY);

            gxxf = magneticXXSpectr + (1 / faceSizeX)^2 .* (1 - forwardShiftX) .* (backwardShiftX - 1) .* electricSpectr;
            gxyf = (1 / (faceSizeX * faceSizeY)) .* (1 - forwardShiftX) .* (backwardShiftY - 1) .* electricSpectr;
            gyxf = (1 / (faceSizeY * faceSizeX)) .* (1 - forwardShiftY) .* (backwardShiftX - 1) .* electricSpectr;
            gyyf = magneticYYSpectr + (1 / faceSizeY)^2 .* (1 - forwardShiftY) .* (backwardShiftY - 1) .* electricSpectr;
        otherwise
            error('Unknown rooftop content: %s', rooftopContent);
    end
end