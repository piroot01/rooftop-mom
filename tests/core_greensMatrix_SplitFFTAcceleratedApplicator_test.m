clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

loggerDir = fullfile(currentDir, '..', '3rd', 'matlab-logger', 'src');
addpath(loggerDir);

% Configure the logger.
logger.Logger.configure( ...
    'timeFormat', 'ss.SSS', ...
    'minLevel', logger.LogLevel.INFO, ...
    'sinks', { ...
        logger.sink.ConsoleSink(), ...
    } ...
);

% Grids covering XONLY, YONLY and BOTH, square and non-square, and sizes with large prime factors.
grids = [20, 1; 1, 23; 6, 6; 4, 11; 13, 5];
sizeX = 1.3;
sizeY = 0.7;

order = 3;
quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);
planner = rooftopMoM.core.greensMatrix.AssemblerChunkPlanner(100 * order^4);
assembler = rooftopMoM.core.greensMatrix.FullAssembler(planner);

k = 0.4;
modifiers = {[], @real, @imag};

for n = 1:size(grids, 1)
    rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(grids(n, 1), grids(n, 2), sizeX, sizeY);

    rooftops = rooftopMoM.basisFunctions.UniformRectPrimitiveUniformRooftops(rectPrimitive, quad);
    rooftops.generate();

    context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);
    greensMatrix = assembler.generate(context);

    x = randn(rooftops.count(), 1) + 1i * randn(rooftops.count(), 1);

    for m = 1:numel(modifiers)
        options = struct('faceSizeX', sizeX, 'faceSizeY', sizeY);
        if ~isempty(modifiers{m})
            options.modifier = modifiers{m};
        end

        fftApplicator = rooftopMoM.core.greensMatrix.FFTAcceleratedApplicator( ...
            context, greensMatrix.generateSpectralGreensMatrix(options), greensMatrix.rooftopContent());
        fftApplicator.generate();

        splitApplicator = rooftopMoM.core.greensMatrix.SplitFFTAcceleratedApplicator( ...
            context, greensMatrix.generateSplitSpectralGreensMatrix(options), greensMatrix.rooftopContent());
        splitApplicator.generate();

        yFFT = fftApplicator.apply(x);
        ySplit = splitApplicator.apply(x);

        relativeError = norm(ySplit - yFFT) / norm(yFFT);
        fprintf('%dx%d, modifier %d: relative error = %.2e\n', grids(n, 1), grids(n, 2), m, relativeError);

        assert(relativeError < 1e-12, 'Split-FFT product differs from the FFT product.');
    end
end
