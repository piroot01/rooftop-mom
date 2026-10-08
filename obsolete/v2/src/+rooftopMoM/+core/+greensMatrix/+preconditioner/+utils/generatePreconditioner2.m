function spectralMatrixData = generatePreconditioner2(generatorHandle, matrixData, faceSizeX, faceSizeY, options)
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

    spatialMatrixData = rooftopMoM.core.greensMatrix.preconditioner.utils.generateSpatialGreensMatrix( ...
        magneticXX, magneticYY, electric, ...
        matrixData.rooftopContent(), faceSizeX, faceSizeY);

    blockNames = {'blockXX', 'blockXY', 'blockYX', 'blockYY'};

    blocks = cell(1, 4);
    for n = 1:numel(blockNames)
        currentBlock = spatialMatrixData.(blockNames{n});

        blockApprox = generatorHandle(currentBlock);

        blocks{n} = fft2(blockApprox);
    end

    spectralMatrixData = rooftopMoM.core.greensMatrix.SpectralMatrixData( ...
        'blockXX', blocks{1}, ...
       'blockXY', blocks{2}, ...
        'blockYX', blocks{3}, ...
        'blockYY', blocks{4} ...
    );
end
