function [gxxf, gxyf, gyxf, gyyf, determinant] = generateSpectralGreensMatrix( ...
        magneticXXCirc, magneticYYCirc, electricCirc, rooftopContent, faceSizeX, faceSizeY)
    switch rooftopContent
        case rooftopMoM.basisFunctions.RooftopContent.XONLY
            magneticXXCircf = fft2(magneticXXCirc);
            clear magneticXXCirc;

            [circSizeX, circSizeY] = size(electricCirc);
            electricCircf = fft2(electricCirc);
            clear electricCirc;

            [forwardShiftX, backwardShiftX, ~, ~] = ...
                rooftopMoM.core.greensMatrix.utils.generateSpectralShiftOperators(circSizeX, circSizeY);

            gxxf = magneticXXCircf + (1 / faceSizeX)^2 .* (1 - forwardShiftX) .* (backwardShiftX - 1) .* electricCircf;
            gxyf = []; gyxf = []; gyyf = [];

            %NOTE: Determinant is the same as gxxf so we wont store its copy.
            determinant = [];
        case rooftopMoM.basisFunctions.RooftopContent.YONLY
            magneticYYCircf = fft2(magneticYYCirc);
            clear magneticYYCirc;

            [circSizeX, circSizeY] = size(electricCirc);
            electricCircf = fft2(electricCirc);
            clear electricCirc;

            [~, ~, forwardShiftY, backwardShiftY] = ...
                rooftopMoM.core.greensMatrix.utils.generateSpectralShiftOperators(circSizeX, circSizeY);

            gyyf = magneticYYCircf + (1 / faceSizeY)^2 .* (1 - forwardShiftY) .* (backwardShiftY - 1) .* electricCircf;
            gxxf = []; gxyf = []; gyxf = [];

            %NOTE: Determinant is the same as gxxf so we wont store its copy.
            determinant = [];
        case rooftopMoM.basisFunctions.RooftopContent.BOTH
            magneticXXCircf = fft2(magneticXXCirc);
            clear magneticXXCirc;

            magneticYYCircf = fft2(magneticYYCirc);
            clear magneticYYCirc;

            [circSizeX, circSizeY] = size(electricCirc);
            electricCircf = fft2(electricCirc);
            clear electricCirc;

            [forwardShiftX, backwardShiftX, forwardShiftY, backwardShiftY] = ...
                rooftopMoM.core.greensMatrix.utils.generateSpectralShiftOperators(circSizeX, circSizeY);

            gxxf = magneticXXCircf + (1 / faceSizeX)^2 .* (1 - forwardShiftX) .* (backwardShiftX - 1) .* electricCircf;
            gxyf = (1 / (faceSizeX * faceSizeY)) .* (1 - forwardShiftX) .* (backwardShiftY - 1) .* electricCircf;
            gyxf = (1 / (faceSizeY * faceSizeX)) .* (1 - forwardShiftY) .* (backwardShiftX - 1) .* electricCircf;
            gyyf = magneticYYCircf + (1 / faceSizeY)^2 .* (1 - forwardShiftY) .* (backwardShiftY - 1) .* electricCircf;
            clear magneticXXCircFFT magneticYYCircFFT electricCircFFT;

            determinant = gxxf .* gyyf - gxyf .* gyxf;
        otherwise
            error('Unknown rooftop content: %s', rooftopContent);
    end
end
