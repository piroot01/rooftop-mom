classdef Problem < handle
% Problem: The electric field integral equation of a Mesh at one wavenumber, discretized with rooftops. The
%   generators of its Green's matrix are assembled once, when an Operator or a Preconditioner needs them first.
%       problem = rooftopMoM.Problem(mesh, 2 * pi, 'quadratureOrder', 5);

    properties (Access = private)
        mesh_ rooftopMoM.Mesh;
        context_ rooftopMoM.assembler.AssemblerContext;

        % rooftopMoM.assembler.greensMatrix.Generators, empty until the first use.
        generators_;
    end

    methods
        function obj = Problem(mesh, wavenumber, options)
        % Problem: Creates the problem; nothing is assembled yet.
        %
        %   Inputs:
        %       mesh: rooftopMoM.Mesh.
        %       wavenumber: Wavenumber in rad/m, 2 * pi / wavelength.
        %
        %   Optional name-value inputs:
        %       'quadratureOrder': Number of Gauss-Legendre points per face and direction (default 5).
        %       'memoryBytes': Memory budget of one chunk of the assembly in bytes (default 2^28).

            arguments
                mesh (1, 1) rooftopMoM.Mesh
                wavenumber (1, 1) {mustBeReal, mustBeFinite, mustBePositive}
                options.quadratureOrder (1, 1) {mustBeInteger, mustBePositive} = 5
                options.memoryBytes (1, 1) {mustBeReal, mustBeFinite, mustBePositive} = 2^28
            end

            quadrature = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(double(options.quadratureOrder));
            rooftops = rooftopMoM.basisFunctions.UniformRooftops(mesh.primitive(), quadrature);

            obj.mesh_ = mesh;
            obj.context_ = rooftopMoM.assembler.AssemblerContext(rooftops, double(wavenumber), ...
                double(options.memoryBytes));

            [countX, countY] = mesh.count();
            [sizeX, sizeY] = mesh.faceSize();
            rooftopMoM.logging.Logger.instance().debug(['%d x %d faces of %g x %g m, %d rooftops, wavenumber %g ' ...
                'rad/m, quadrature order %d.'], countX, countY, sizeX, sizeY, rooftops.count(), wavenumber, ...
                options.quadratureOrder);
        end
    end

    methods
        function mesh = mesh(obj)
            mesh = obj.mesh_;
        end

        function wavenumber = wavenumber(obj)
            wavenumber = obj.context_.wavenumber();
        end

        function count = rooftopCount(obj)
            count = obj.context_.rooftops().count();
        end

        function rooftops = rooftops(obj)
            rooftops = obj.context_.rooftops();
        end

        function context = context(obj)
            context = obj.context_;
        end

        function generators = generators(obj)
        % generators: Returns the generators of the Green's matrix, assembled on the first call.

            if isempty(obj.generators_)
                obj.generators_ = rooftopMoM.assembler.greensMatrix.FullAssembler().generate(obj.context_);
            end

            generators = obj.generators_;
        end

        function bytes = memory(obj)
        % memory: Returns the bytes of the generators, zero before they are assembled.

            bytes = 0;

            if ~isempty(obj.generators_)
                bytes = obj.generators_.memory();
            end
        end
    end
end
