function impedanceMatrix = generateImpedanceMatrix(varargin)
% generateImpedanceMatrix: Assembles the full (dense) impedance matrix of a uniform rectangular antenna using the
%   interaction-symmetries accelerated assembler.
%
%   Inputs:
%       varargin: Name-value pairs:
%           'countX': Number of faces in the x direction (default 10).
%           'countY': Number of faces in the y direction (default 1).
%           'sizeX': Face size in the x direction (default 1).
%           'sizeY': Face size in the y direction (default 1).
%           'quadratureOrder': Order of the Gauss-Legendre quadrature (default 5).
%           'wavenumber': Wavenumber of the excitation (default 2 * pi / 10).
%           'chunkSize': Number of the distance matrix elements computed per chunk (default 1e7).
%
%   Output:
%       impedanceMatrix: NxN complex impedance matrix with the rooftops ordered as in prepareSolverContext.

    p = inputParser();

    addParameter(p, 'countX', 10);
    addParameter(p, 'countY', 1);
    addParameter(p, 'sizeX', 1);
    addParameter(p, 'sizeY', 1);
    addParameter(p, 'quadratureOrder', 5);
    addParameter(p, 'wavenumber', 2 * pi / 10);
    addParameter(p, 'chunkSize', 1e7);

    parse(p, varargin{:});

    % Create the mesh of the antenna
    rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(p.Results.countX, p.Results.countY, p.Results.sizeX, p.Results.sizeY);
    mesh = rooftopMoM.mesh.UniformRect(rectPrimitive).generateMesh();

    % Generate rooftop basis functions on the mesh
    quadrature = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(p.Results.quadratureOrder);

    rooftops = rooftopMoM.basisFunctions.MeshUniformRooftops(mesh, quadrature);
    rooftops.generate();

    % Assemble the full impedance matrix reusing the symmetric rooftop interactions
    planner = rooftopMoM.core.impedanceMatrix.AssemblerChunkPlanner(p.Results.chunkSize);
    context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quadrature, rooftops, p.Results.wavenumber);

    fullAssembler = rooftopMoM.core.impedanceMatrix.FullAssembler(planner);
    accelerator = rooftopMoM.core.impedanceMatrix.InteractionSymmetriesAccelerator(fullAssembler);

    impedanceMatrix = accelerator.generateFull(context).value();

    TRACE("generateImpedanceMatrix: Impedance matrix assembled: rooftop count = %g", rooftops.count());
end
