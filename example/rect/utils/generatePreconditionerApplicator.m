function preconditionerApplicator = generatePreconditionerApplicator(solverContext, varargin)
% generatePreconditionerApplicator: Generates the FFT accelerated applicator of a circulant preconditioner of the
%   greens matrix.
%
%   Inputs:
%       solverContext: Struct from prepareSolverContext with the fields context and greensMatrix.
%       varargin: Name-value pairs:
%           'sizeX': Face size in the x direction (default 1).
%           'sizeY': Face size in the y direction (default 1).
%           'preconditionerType': Circulant approximation, 'TChan' or 'Strang' (default 'TChan').
%           'modifier': Optional function handle applied to the greens matrix parts, e.g. @imag (default none).
%
%   Output:
%       preconditionerApplicator: Generated preconditioner applicator.

    p = inputParser();

    addParameter(p, 'sizeX', 1);
    addParameter(p, 'sizeY', 1);
    addParameter(p, 'preconditionerType', 'TChan');
    addParameter(p, 'modifier', []);

    parse(p, varargin{:});

    switch p.Results.preconditionerType
        case 'TChan'
            preconditionerHandle = @rooftopMoM.core.greensMatrix.preconditioner.utils.generateTChanPreconditionerKernel;
        case 'Strang'
            preconditionerHandle = @rooftopMoM.core.greensMatrix.preconditioner.utils.generateStrangPreconditionerKernel;
        otherwise
            error('Unknown preconditioner type: %s', p.Results.preconditionerType);
    end

    options = struct();

    if ~isempty(p.Results.modifier)
        options.modifier = p.Results.modifier;
    end

    preconditioner = rooftopMoM.core.greensMatrix.preconditioner.utils.generatePreconditioner1( ...
        preconditionerHandle, solverContext.greensMatrix, p.Results.sizeX, p.Results.sizeY, options);

    preconditionerApplicator = rooftopMoM.core.greensMatrix.preconditioner.FFTAcceleratedApplicator( ...
        solverContext.context, preconditioner, solverContext.greensMatrix.rooftopContent());
    preconditionerApplicator.generate();
end
