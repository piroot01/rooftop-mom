function spectralMatrixData = generatePreconditioner3(generatorHandle, matrixData, faceSizeX, faceSizeY, options)
    if nargin < 3
        faceSizeX = 1;
        faceSizeY = 1;
    end

    if nargin > 4 && isfield(options, 'modifier')
        magneticXX = options.modifier(matrixData.magneticXX());
        magneticYY = options.modifier(matrixData.magneticYY());
        electric = options.modifier(matrixData.electric());
    else
        magneticXX = matrixData.magneticXX();
        magneticYY = matrixData.magneticYY();
        electric = matrixData.electric();
    end

    % Same direction handled the same way as in the generatePreconditioner2
    [blockXX, blockYY] = rooftopMoM.core.greensMatrix.preconditioner.utils.generateSameDirSpatialGreensMatrix( ...
        magneticXX, magneticYY, electric, ...
        matrixData.rooftopContent(), faceSizeX, faceSizeY);

    blockXX = fft2(generatorHandle(blockXX));
    blockYY = fft2(generatorHandle(blockYY));

    % Cross direction interactions handled the same way as in the generatePreconditioner1
    electricApprox = generatorHandle(electric);

    [blockXY, blockYX] = rooftopMoM.core.greensMatrix.preconditioner.utils.generateCrossDirSpectralGreensMatrix( ...
        electricApprox, ...
        faceSizeX, faceSizeY);

    spectralMatrixData = rooftopMoM.core.greensMatrix.SpectralMatrixData( ...
        'blockXX', blockXX, ...
        'blockXY', blockXY, ...
        'blockYX', blockYX, ...
        'blockYY', blockYY ...
    );
end
