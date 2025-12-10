function output = generateApplicators(solverContext, varargin)
    p = inputParser();

    addParameter(p, 'sizeX', 1);
    addParameter(p, 'sizeY', 1);
    addParameter(p, 'preconditionerType', 'TChan');

    parse(p, varargin{:});

    sizeX = p.Results.sizeX;
    sizeY = p.Results.sizeY;
    preconditionerType = p.Results.preconditionerType;

    % real(Z) applicator
    spectralRealGreensMatrix = solverContext.greensMatrix.generateSpectralGreensMatrix(struct('modifier', @real));

    realGreensMatrixApplicator = rooftopMoM.core.greensMatrix.FFTAcceleratedApplicator( ...
        solverContext.context, spectralRealGreensMatrix, solverContext.greensMatrix.rooftopContent());

    % imag(Z) applicator
    spectralImagGreensMatrix = solverContext.greensMatrix.generateSpectralGreensMatrix(struct('modifier', @imag));

    imagGreensMatrixApplicator = rooftopMoM.core.greensMatrix.FFTAcceleratedApplicator( ...
        solverContext.context, spectralImagGreensMatrix, solverContext.greensMatrix.rooftopContent());

    if strcmp(preconditionerType, 'TChan')
        preconditionerHandle = @rooftopMoM.core.greensMatrix.preconditioner.utils.generateTChanPreconditionerKernel;
    end

    preconditioner = rooftopMoM.core.greensMatrix.preconditioner.utils.generatePreconditioner1( ...
        preconditionerHandle, ...
        solverContext.greensMatrix, sizeX, sizeY, ...
        struct('modifier', @imag));

    preconditionerApplicator = rooftopMoM.core.greensMatrix.preconditioner.FFTAcceleratedApplicator( ...
        solverContext.context, preconditioner, solverContext.greensMatrix.rooftopContent());
    preconditionerApplicator.generate();

    output = struct( ...
        'realGreensMatrix', realGreensMatrixApplicator, ...
        'imagGreensMatrix', imagGreensMatrixApplicator, ...
        'imagPreconditioner', preconditionerApplicator ...
    );
end
