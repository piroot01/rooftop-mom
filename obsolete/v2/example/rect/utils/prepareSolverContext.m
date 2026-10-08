function context = prepareSolverContext(varargin)
    p = inputParser();

    addParameter(p, 'countX', 10);
    addParameter(p, 'countY', 1);
    addParameter(p, 'sizeX', 1);
    addParameter(p, 'sizeY', 1);
    addParameter(p, 'quadratureOrder', 5);
    addParameter(p, 'wavenumber', 2 * pi / 10);
    addParameter(p, 'chunkSize', 1e7);
    addParameter(p, 'waveDirection', [0; 0; -1]);
    addParameter(p, 'wavePolarization', [0; 1; 0]);
    addParameter(p, 'waveIntensity', 1);

    parse(p, varargin{:});

    countX = p.Results.countX;
    countY = p.Results.countY;
    sizeX = p.Results.sizeX;
    sizeY = p.Results.sizeY;
    quadratureOrder = p.Results.quadratureOrder;
    wavenumber = p.Results.wavenumber;
    chunkSize = p.Results.chunkSize;
    waveDirection = p.Results.waveDirection;
    wavePolarization = p.Results.wavePolarization;
    waveIntensity = p.Results.waveIntensity;

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

    TRACE("solveAntennaCurrent: Rooftop basis functions generated: rooftop count = %g", rooftops.count());

    % Prepare a Green's matrix filling strategy and assemble the context
    planner = rooftopMoM.core.greensMatrix.AssemblerChunkPlanner(chunkSize);
    context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quadrature, rooftops, wavenumber);

    % Instantiate a full Green's matrix assembler 
    assembler = rooftopMoM.core.greensMatrix.FullAssembler(planner);

    % Generate the Green's matrix
    greensMatrix = assembler.generate(context);

    TRACE("solveAntennaCurrent: Green's matrix assembled.");

    % Generate the excitation vector using the cached rooftop data.
    %excitation = rooftopMoM.operator.PlaneWaveExcitation.generateWithAssemblerContext(context, waveDirection, wavePolarization, waveIntensity, chunkSize);
    excitation = rooftopMoM.operator.PlaneWaveExcitation.generatePerpendicular(rooftops, wavePolarization, waveIntensity);

    TRACE("solveAntennaCurrent: Plane wave excitation vector assembled.");

    context = struct( ...
        'context', context, ...
        'greensMatrix', greensMatrix, ...
        'excitation', excitation ...
    );
end
