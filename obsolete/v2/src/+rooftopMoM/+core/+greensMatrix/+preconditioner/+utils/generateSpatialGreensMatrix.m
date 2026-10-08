function spectralMatrixData = generateSpatialGreensMatrix(magneticXX, magneticYY, electric, rooftopContent, faceSizeX, faceSizeY)
    if ~isempty(electric)
        [kernelSizeX, kernelSizeY] = size(electric);
    elseif ~isempty(magneticXX)
        [kernelSizeX, kernelSizeY] = size(magneticXX);
    elseif ~isempty(magneticYY)
        [kernelSizeX, kernelSizeY] = size(magneticYY);
    else
        error('At least one of the kernels must be nonempty!');
    end

    isXActive = (rooftopContent == rooftopMoM.basisFunctions.RooftopContent.XONLY) || ...
        (rooftopContent == rooftopMoM.basisFunctions.RooftopContent.BOTH);

    isYActive = (rooftopContent == rooftopMoM.basisFunctions.RooftopContent.YONLY) || ...
        (rooftopContent == rooftopMoM.basisFunctions.RooftopContent.BOTH);

    blockXX = []; blockXY = []; blockYX = []; blockYY = [];

    if isempty(magneticXX)
        magneticXX = 0;
    end

    if isempty(magneticYY)
        magneticYY = 0;
    end

    if isempty(electric)
        electric = 0;
    end

    if isXActive
        termElectric = zeros([kernelSizeX, kernelSizeY], 'like', electric);

        termElectric(2:(end - 1), :) = 2 * electric(2:(end - 1), :) - electric(1:(end - 2), :) - electric(3:end, :);
        termElectric(1, :) = 2 * (electric(1, :) - electric(2, :));
        termElectric(end, :) = 2 * (electric(end, :) - electric(end - 1, :));

        termMagnetic = rooftopMoM.core.greensMatrix.preconditioner.utils.padWithZeros(magneticXX, kernelSizeX, kernelSizeY);

        blockXX = termMagnetic - (1 / faceSizeX)^2 .* termElectric;
    end

    if isYActive
        termElectric = zeros([kernelSizeX, kernelSizeY], 'like', electric);

        termElectric(:, 2:(end - 1)) = 2 * electric(:, 2:(end - 1)) - electric(:, 1:(end - 2)) - electric(:, 3:end);
        termElectric(:, 1) = 2 * (electric(:, 1) - electric(:, 2));
        termElectric(:, end) = 2 * (electric(:, end) - electric(:, end - 1));

        termMagnetic = rooftopMoM.core.greensMatrix.preconditioner.utils.padWithZeros(magneticYY, kernelSizeX, kernelSizeY);

        blockYY = termMagnetic - (1 / faceSizeY)^2 .* termElectric;
    end

    if isXActive && isYActive
        % Corresponds to: T(i,j) = E(i, j-1) - E(i, j)
        diffY = zeros(size(electric), 'like', electric);
        diffY(:, 2:end) = electric(:, 1:(end - 1)) - electric(:, 2:end);
        diffY(:, 1) = (-1) * diffY(:, 2);

        % Corresponds to: R(i,j) = T(i,j) - T(i+1, j)
        termElectric = zeros(size(electric), 'like', electric);
        termElectric(1:(end - 1), :) = diffY(1:(end - 1), :) - diffY(2:end, :);
        termElectric(end, :) = (-1) * termElectric(end - 1, :);

        blockXY = (1 / (faceSizeX * faceSizeY)) .* termElectric;

        % Corresponds to: T(i,j) = E(i-1, j) - E(i, j)
        diffX = zeros(size(electric), 'like', electric);
        diffX(2:end, :) = electric(1:(end - 1), :) - electric(2:end, :);
        diffX(1, :) = (-1) * diffX(2, :);

        % Corresponds to: R(i,j) = T(i,j) - T(i, j+1)
        termElectric = zeros(size(electric), 'like', electric);
        termElectric(:, 1:(end - 1)) = diffX(:, 1:(end - 1)) - diffX(:, 2:end);
        termElectric(:, end) = (-1) * termElectric(:, end - 1);

        blockYX = (1 / (faceSizeY * faceSizeX)) .* termElectric;
    end

    spectralMatrixData = rooftopMoM.core.greensMatrix.SpectralMatrixData( ...
        'blockXX', blockXX, ...
        'blockXY', blockXY, ...
        'blockYX', blockYX, ...
        'blockYY', blockYY ...
    );
end
