function spectralMatrixData = generatePreconditioner1(generatorHandle, matrixData, faceSizeX, faceSizeY, options)
    if nargin < 3
        faceSizeX = 1;
        faceSizeY = 1;
    end

    fieldNames = {'magneticXX', 'magneticYY', 'electric'};

    circulantApprox = cell(1, 3);
    for n = 1:numel(fieldNames)
        currentPart = matrixData.(fieldNames{n});

        if nargin > 4 && isfield(options, 'modifier')
            currentPart = options.modifier(currentPart);
        end

        [kernelSizeX, kernelSizeY] = matrixData.minimalKernelSize();

        currentPartPadded = rooftopMoM.core.greensMatrix.preconditioner.utils.padWithZeros(currentPart, kernelSizeX, kernelSizeY);

        circulantApprox{n} = generatorHandle(currentPartPadded);
    end

    spectralMatrixData = rooftopMoM.core.greensMatrix.utils.generateSpectralGreensMatrix( ...
        circulantApprox{1}, circulantApprox{2}, circulantApprox{3}, matrixData.rooftopContent(), faceSizeX, faceSizeY);
end
