function [blockXY, blockYX] = generateCrossDirSpectralGreensMatrix(electricCirc, faceSizeX, faceSizeY)
    electricCircSpectr = fft2(electricCirc);

    if isempty(electricCircSpectr)
        blockXY = []; blockYX = [];
        return;
    end

    [circSizeX, circSizeY] = size(electricCirc);

    % Generate Shift Operators (Calculate once, use everywhere)
    [forwardShiftX, backwardShiftX, forwardShiftY, backwardShiftY] = ...
        rooftopMoM.core.greensMatrix.utils.generateSpectralShiftOperators(circSizeX, circSizeY);

    blockXY = (1 / (faceSizeX * faceSizeY)) .* (1 - forwardShiftX) .* (backwardShiftY - 1) .* electricCircSpectr;
    blockYX = (1 / (faceSizeY * faceSizeX)) .* (1 - forwardShiftY) .* (backwardShiftX - 1) .* electricCircSpectr;
end
