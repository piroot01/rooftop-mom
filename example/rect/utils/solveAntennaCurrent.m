function output = solveAntennaCurrent(varargin)
    p = inputParser();

    addParameter(p, 'countX', 10);
    addParameter(p, 'countY', 1);
    addParameter(p, 'sizeX', 1);
    addParameter(p, 'sizeY', 1);
    addParameter(p, 'quadratureOrder', 5);
    addParameter(p, 'wavenumber', 2 * pi / 10);
    addParameter(p, 'elementCount', 1e7);
    addParameter(p, 'excitation', [0; 0; 0; 0; 1; 0; 0; 0; 0]);
    addParameter(p, 'tolerance', 1e-12);
    addParameter(p, 'maxIterCount', 1000);
    addParameter(p, 'useSplitFFT', false);
    addParameter(p, 'genome', []);

    parse(p, varargin{:});

    countX = p.Results.countX;
    countY = p.Results.countY;
    sizeX = p.Results.sizeX;
    sizeY = p.Results.sizeY;
    quadratureOrder = p.Results.quadratureOrder;
    wavenumber = p.Results.wavenumber;
    elementCount = p.Results.elementCount;
    excitation = p.Results.excitation;
    tolerance = p.Results.tolerance;
    maxIterCount = p.Results.maxIterCount;
    useSplitFFT = p.Results.useSplitFFT;

    DEBUG("solveAntennaCurrent: Parsed inputs: countX = %g, countY = %g, sizeX = %g m, sizeY = %g m.", countX, countY, sizeX, sizeY);

    % create the mesh of the antenna
    rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);

    TRACE("solveAntennaCurrent: Rectangle primitive generated.");

    % Prepare the quadrature points and weights
    quadrature = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(quadratureOrder);

    TRACE("solveAntennaCurrent: Quadrature points and weights generated.");

    % Generate rooftop basis functions
    rooftops = rooftopMoM.basisFunctions.UniformRectPrimitiveUniformRooftops(rectPrimitive, quadrature);
    rooftops.generate();

    TRACE("solveAntennaCurrent: Rooftop basis functions generated.");

    % Prepare a Green's matrix filling strategy and assemble the context
    planner = rooftopMoM.core.greensMatrix.AssemblerChunkPlanner(elementCount);
    context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quadrature, rooftops, wavenumber);

    % Instantiate a full Green's matrix assembler 
    assembler = rooftopMoM.core.greensMatrix.FullAssembler(planner);

    % Generate the Green's matrix
    greensMatrix = assembler.generateFull(context);

    TRACE("solveAntennaCurrent: Green's matrix assembled.");

    % Create the Green's matrix applicator
    if ~useSplitFFT
        greensMatrixApplicator = rooftopMoM.core.greensMatrix.FFTAcceleratedApplicator(context, greensMatrix);
    else
        greensMatrixApplicator = rooftopMoM.core.greensMatrix.SplitFFTAcceleratedApplicator(context, greensMatrix);
    end

    greensMatrixApplicator.generate();

    TRACE("solveAntennaCurrent: Green's matrix applicator generated.");

    % Instantiate the preconditioned conjugate gradient solver
    solver = rooftopMoM.core.solver.SpectrallyPreconditionedConjugateGradient(greensMatrixApplicator);

    % Finally solve for the current
    solverOutput = solver.solve( ...
        'excitation', excitation, ...
        'tolerance', tolerance, ...
        'maxIterCount', maxIterCount ...
    );

    output.solverOutput = solverOutput;
end
