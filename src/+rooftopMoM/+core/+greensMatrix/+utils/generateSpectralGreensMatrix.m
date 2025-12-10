function spectralMatrixData = generateSpectralGreensMatrix(magneticXXCirc, magneticYYCirc, electricCirc, rooftopContent, faceSizeX, faceSizeY)
    if ~isempty(electricCirc)
        [circSizeX, circSizeY] = size(electricCirc);
    elseif ~isempty(magneticXXCirc)
        [circSizeX, circSizeY] = size(magneticXXCirc);
    elseif ~isempty(magneticYYCirc)
        [circSizeX, circSizeY] = size(magneticYYCirc);
    else
        error('At least one of the circulants must be nonempty!');
    end

    % Compute Spectra (FFT) & Sanitize
    magneticXXCircSpectr = fft2(magneticXXCirc);
    if isempty(magneticXXCircSpectr)
        magneticXXCircSpectr = 0;
    end

    magneticYYCircSpectr = fft2(magneticYYCirc);
    if isempty(magneticYYCircSpectr)
        magneticYYCircSpectr = 0;
    end

    electricCircSpectr = fft2(electricCirc);
    if isempty(electricCircSpectr)
        electricCircSpectr = 0;
    end

    % Generate Shift Operators (Calculate once, use everywhere)
    [forwardShiftX, backwardShiftX, forwardShiftY, backwardShiftY] = ...
        rooftopMoM.core.greensMatrix.utils.generateSpectralShiftOperators(circSizeX, circSizeY);

    % Determine Active Modes
    isXActive = (rooftopContent == rooftopMoM.basisFunctions.RooftopContent.XONLY) || ...
        (rooftopContent == rooftopMoM.basisFunctions.RooftopContent.BOTH);

    isYActive = (rooftopContent == rooftopMoM.basisFunctions.RooftopContent.YONLY) || ...
        (rooftopContent == rooftopMoM.basisFunctions.RooftopContent.BOTH);

    % Compute Blocks
    blockXX = []; blockXY = []; blockYX = []; blockYY = [];

    % Compute XX Component
    if isXActive
        termMagnetic = magneticXXCircSpectr;
        termElectric = (1 / faceSizeX)^2 .* (1 - forwardShiftX) .* (backwardShiftX - 1) .* electricCircSpectr;
        blockXX = termMagnetic + termElectric;
    end

    % Compute YY Component
    if isYActive
        termMagnetic = magneticYYCircSpectr;
        termElectric = (1 / faceSizeY)^2 .* (1 - forwardShiftY) .* (backwardShiftY - 1) .* electricCircSpectr;
        blockYY = termMagnetic + termElectric;
    end

    % Compute Cross-Coupling (XY / YX)
    if isXActive && isYActive
        blockXY = (1 / (faceSizeX * faceSizeY)) .* (1 - forwardShiftX) .* (backwardShiftY - 1) .* electricCircSpectr;
        blockYX = (1 / (faceSizeY * faceSizeX)) .* (1 - forwardShiftY) .* (backwardShiftX - 1) .* electricCircSpectr;
    end

    spectralMatrixData = rooftopMoM.core.greensMatrix.SpectralMatrixData( ...
        'blockXX', blockXX, ...
        'blockXY', blockXY, ...
        'blockYX', blockYX, ...
        'blockYY', blockYY ...
    );
end
