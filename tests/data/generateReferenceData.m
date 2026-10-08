function generateReferenceData()
% generateReferenceData: Computes the reference data of the v2 code in obsolete/v2 and saves it as reference.mat next
%   to this file. The reference tests of v3 (tests/+reference) compare against it.
%   Run it in a fresh MATLAB session, so the v3 package is not on the path:
%       matlab -batch "cd tests/data; generateReferenceData"

    here = fileparts(mfilename('fullpath'));
    root = fileparts(fileparts(here));
    v2 = fullfile(root, 'obsolete', 'v2');

    % Both versions are the package rooftopMoM, so the v3 source must not be on the path.
    if contains(path, fullfile(root, 'src'))
        rmpath(fullfile(root, 'src'));
    end
    addpath(fullfile(v2, 'src'), fullfile(v2, '3rd', 'matlab-logger', 'src'));
    logger.Logger.configure('minLevel', logger.LogLevel.OFF);

    rng(1);

    % countX, countY, sizeX, sizeY, quadrature order, wavenumber
    cases = {
        5, 4, 0.1, 0.15, 3, 2 * pi
        6, 1, 0.1, 0.02, 4, 2 * pi
        1, 5, 0.03, 0.1, 3, 2 * pi
    };

    for n = size(cases, 1):(-1):1
        reference(n) = generateCase_(cases{n, :});
    end

    reference(1).impedanceMatrix = generateImpedanceMatrix_(cases{1, :});

    % The dense matrix comes from a separate code path of v2; it should agree with the fast product.
    difference = norm(reference(1).impedanceMatrix * reference(1).vector - reference(1).productFFT) / ...
        norm(reference(1).productFFT);
    fprintf('v2 dense against FFT product: relative difference %.1e\n', difference);

    save(fullfile(here, 'reference.mat'), 'reference');
end

function data = generateCase_(countX, countY, sizeX, sizeY, order, wavenumber)
    primitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
    quadrature = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);
    rooftops = rooftopMoM.basisFunctions.UniformRectPrimitiveUniformRooftops(primitive, quadrature);
    rooftops.generate();
    context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quadrature, rooftops, wavenumber);
    planner = rooftopMoM.core.greensMatrix.AssemblerChunkPlanner(1e6);
    matrixData = rooftopMoM.core.greensMatrix.FullAssembler(planner).generate(context);
    content = matrixData.rooftopContent();
    options = struct('faceSizeX', sizeX, 'faceSizeY', sizeY);
    vector = complex(randn(rooftops.count(), 1), randn(rooftops.count(), 1));

    data.countX = countX;
    data.countY = countY;
    data.sizeX = sizeX;
    data.sizeY = sizeY;
    data.order = order;
    data.wavenumber = wavenumber;
    data.magneticXX = matrixData.magneticXX();
    data.magneticYY = matrixData.magneticYY();
    data.electric = matrixData.electric();
    data.vector = vector;

    data.productFFT = apply_(rooftopMoM.core.greensMatrix.FFTAcceleratedApplicator(context, ...
        matrixData.generateSpectralGreensMatrix(options), content), vector);
    data.productSplit = apply_(rooftopMoM.core.greensMatrix.SplitFFTAcceleratedApplicator(context, ...
        matrixData.generateSplitSpectralGreensMatrix(options), content), vector);
    data.productLean = apply_(rooftopMoM.core.greensMatrix.MemoryEfficientSplitFFTAcceleratedApplicator(context, ...
        matrixData.generateSplitSpectralGreensMatrix(options), content), vector);

    % Rows: the Strang, T. Chan and R. Chan kernel; columns: generatePreconditioner1, 2 and 3.
    kernels = {
        @rooftopMoM.core.greensMatrix.preconditioner.utils.generateStrangPreconditionerKernel
        @rooftopMoM.core.greensMatrix.preconditioner.utils.generateTChanPreconditionerKernel
        @rooftopMoM.core.greensMatrix.preconditioner.utils.generateRChanPreconditionerKernel
    };
    builders = {
        @rooftopMoM.core.greensMatrix.preconditioner.utils.generatePreconditioner1
        @rooftopMoM.core.greensMatrix.preconditioner.utils.generatePreconditioner2
        @rooftopMoM.core.greensMatrix.preconditioner.utils.generatePreconditioner3
    };
    data.preconditioned = cell(numel(kernels), numel(builders));

    for k = 1:numel(kernels)
        for b = 1:numel(builders)
            spectra = builders{b}(kernels{k}, matrixData, sizeX, sizeY);
            preconditioner = rooftopMoM.core.greensMatrix.preconditioner.FFTAcceleratedApplicator(context, ...
                spectra, content);
            data.preconditioned{k, b} = preconditioner.apply(vector);
        end
    end

    % An oblique plane wave: the direction and the polarization are perpendicular.
    data.direction = [1; 2; -2] / 3;
    data.polarization = [2; -1; 0] / sqrt(5);
    data.intensity = 2;
    data.planeWave = rooftopMoM.operator.PlaneWaveExcitation.generate(quadrature, rooftops, wavenumber, ...
        data.direction, data.polarization, data.intensity, 1e6);

    data.impedanceMatrix = [];
end

function output = apply_(applicator, vector)
    applicator.generate();
    output = applicator.apply(vector);
end

function impedanceMatrix = generateImpedanceMatrix_(countX, countY, sizeX, sizeY, order, wavenumber)
    primitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
    mesh = rooftopMoM.mesh.UniformRect(primitive).generateMesh();
    quadrature = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);
    rooftops = rooftopMoM.basisFunctions.MeshUniformRooftops(mesh, quadrature);
    rooftops.generate();
    context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quadrature, rooftops, wavenumber);
    planner = rooftopMoM.core.impedanceMatrix.AssemblerChunkPlanner(1e6);
    % The plain dense assembler of v2 fails for non-square faces; the symmetry accelerator handles them.
    fullAssembler = rooftopMoM.core.impedanceMatrix.FullAssembler(planner);
    impedanceMatrix = rooftopMoM.core.impedanceMatrix.InteractionSymmetriesAccelerator(fullAssembler) ...
        .generateFull(context).value();
end
