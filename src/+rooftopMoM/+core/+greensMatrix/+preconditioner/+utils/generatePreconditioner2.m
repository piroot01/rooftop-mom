function spectralMatrixData = generatePreconditioner2(generatorHandle, matrixData, faceSizeX, faceSizeY)
    if nargin < 3
        faceSizeX = 1;
        faceSizeY = 1;
    end

    electric = matrixData.electric();
    blockXX = zeros(size(electric), 'like', electric);

    blockXX(2:(end - 1), :) = 2 * electric(2:(end - 1), :) - electric(1:(end - 2), :) - electric(3:end, :);
    blockXX(1, :) = 2 * (electric(1, :) - electric(2, :));
    blockXX(end, :) = 2 * (electric(end, :) - electric(end - 1, :));

    blockXXApprox = generatorHandle(blockXX);

    spectralMatrixData = rooftopMoM.core.greensMatrix.SpectralMatrixData( ...
        'blockXX', blockXXApprox, ...
        'blockXY', [], ...
        'blockYX', [], ...
        'blockYY', [] ...
    );
end
