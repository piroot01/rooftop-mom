function greensMatrixApplicator = generateGreensMatrixApplicator(solverContext, varargin)
% generateGreensMatrixApplicator: Generates the FFT or the split-FFT accelerated applicator of the greens matrix.
%
%   Inputs:
%       solverContext: Struct from prepareSolverContext with the fields context and greensMatrix.
%       varargin: Name-value pairs:
%           'sizeX': Face size in the x direction (default 1).
%           'sizeY': Face size in the y direction (default 1).
%           'useSplitFFT': Logical flag selecting the split-FFT applicator instead of the FFT one (default false).
%           'memoryEfficient': Logical flag selecting the memory efficient split-FFT applicator (default false).
%           'modifier': Optional function handle applied to the greens matrix parts, e.g. @real (default none).
%
%   Output:
%       greensMatrixApplicator: Generated greens matrix applicator.

    p = inputParser();

    addParameter(p, 'sizeX', 1);
    addParameter(p, 'sizeY', 1);
    addParameter(p, 'useSplitFFT', false);
    addParameter(p, 'memoryEfficient', false);
    addParameter(p, 'modifier', []);

    parse(p, varargin{:});

    options = struct('faceSizeX', p.Results.sizeX, 'faceSizeY', p.Results.sizeY);

    if ~isempty(p.Results.modifier)
        options.modifier = p.Results.modifier;
    end

    greensMatrix = solverContext.greensMatrix;

    if p.Results.useSplitFFT && p.Results.memoryEfficient
        greensMatrixApplicator = rooftopMoM.core.greensMatrix.MemoryEfficientSplitFFTAcceleratedApplicator( ...
            solverContext.context, greensMatrix.generateSplitSpectralGreensMatrix(options), greensMatrix.rooftopContent());
    elseif p.Results.useSplitFFT
        greensMatrixApplicator = rooftopMoM.core.greensMatrix.SplitFFTAcceleratedApplicator( ...
            solverContext.context, greensMatrix.generateSplitSpectralGreensMatrix(options), greensMatrix.rooftopContent());
    else
        greensMatrixApplicator = rooftopMoM.core.greensMatrix.FFTAcceleratedApplicator( ...
            solverContext.context, greensMatrix.generateSpectralGreensMatrix(options), greensMatrix.rooftopContent());
    end

    greensMatrixApplicator.generate();
end
