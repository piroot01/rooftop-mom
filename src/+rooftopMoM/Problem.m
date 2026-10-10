classdef Problem < handle
% Problem: The electric field integral equation of a Mesh at one wavenumber, discretized with rooftops. The
%   generators of its Green's matrix are assembled once, when an Operator or a Preconditioner needs them first.
%       problem = rooftopMoM.Problem(mesh, 2 * pi, 'quadratureOrder', [16, 8]);
%       problem = rooftopMoM.Problem(mesh, 'electricalSize', 4);
%       problem = rooftopMoM.Problem.fromElectricalSize(4, 'density', 10, 'aspect', 2);
%   The current depends only on the electrical size ka, with a the radius of the smallest sphere around the plate, and
%   on the shape and the faces per wavelength of the mesh: a plate scaled together with the wavelength gives the same
%   current, input impedance and radar cross section over the squared wavelength.

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
        %       wavenumber: Wavenumber in rad/m, 2 * pi / wavelength; left out when 'electricalSize' is given.
        %
        %   Optional name-value inputs:
        %       'electricalSize': ka in place of the wavenumber, which becomes ka / mesh.radius().
        %       'quadratureOrder': Gauss-Legendre orders [near, far] of the cell pairs, or one order for both (default
        %           [16, 8], which gives the entries of the Green's matrix to about 1e-14, or 1e-13 for faces ten times
        %           as long as wide). The faces of a near pair are less than three longer face sides apart. The
        %           excitation uses the far order.
        %       'chunkSize': Bytes of one chunk of the assembly, its peak memory (default the chunk size set in
        %           rooftopMoM.Resources.configure, otherwise 2^28 bytes or what the memory ceiling leaves), at most
        %           what the assembly can use: its near or its far pass in one chunk. Every chunk costs a fixed time,
        %           so chunks below 2^24 bytes are slower (2^20 takes about 1.7 times as long); larger ones are as fast
        %           or faster.

            arguments
                mesh (1, 1) rooftopMoM.Mesh
                wavenumber {mustBeScalarOrEmpty, mustBeReal, mustBeFinite, mustBePositive} = []
                options.electricalSize {mustBeScalarOrEmpty, mustBeReal, mustBeFinite, mustBePositive} = []
                options.quadratureOrder {mustBeInteger, mustBePositive, mustBeOneOrTwo_} = [16, 8]
                options.chunkSize {mustBeScalarOrEmpty, mustBeReal, mustBeFinite, mustBePositive} = []
            end

            if isempty(wavenumber) == isempty(options.electricalSize)
                error('rooftopMoM:problem:invalidWavenumber', ...
                    'A problem takes either the wavenumber or the ''electricalSize'', not both or neither.');
            end

            if isempty(wavenumber)
                wavenumber = options.electricalSize / mesh.radius();
            end

            orders = double([options.quadratureOrder(1), options.quadratureOrder(end)]);
            quadrature = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(orders(2));
            rooftops = rooftopMoM.basisFunctions.UniformRooftops(mesh.primitive(), quadrature);
            chunkSize = double(options.chunkSize);

            if isempty(chunkSize)
                chunkSize = rooftopMoM.resources.Budget.instance().chooseChunk( ...
                    rooftopMoM.resources.estimateGenerators(rooftops));
            end

            % A larger chunk would change nothing, so the memory checks count only what the assembly can use.
            chunkSize = min(chunkSize, rooftopMoM.resources.estimateChunk(rooftops, orders));

            obj.mesh_ = mesh;
            obj.context_ = rooftopMoM.assembler.AssemblerContext(rooftops, double(wavenumber), chunkSize, orders);
            obj.requireAssembly_();

            [countX, countY] = mesh.count();
            [sizeX, sizeY] = mesh.faceSize();
            rooftopMoM.logging.Logger.instance().debug(['%d x %d faces of %g x %g m, %d rooftops, wavenumber %g ' ...
                'rad/m, ka %.4g, %.3g faces per wavelength, quadrature orders %d near and %d far.'], countX, countY, ...
                sizeX, sizeY, rooftops.count(), wavenumber, obj.electricalSize(), obj.density(), orders(1), orders(2));
        end

        function delete(obj)
        % delete: Returns the memory of the generators to the budget of rooftopMoM.Resources.

            rooftopMoM.resources.Budget.instance().release(obj.memory());
        end
    end

    methods
        function mesh = mesh(obj)
            mesh = obj.mesh_;
        end

        function wavenumber = wavenumber(obj)
            wavenumber = obj.context_.wavenumber();
        end

        function electricalSize = electricalSize(obj)
        % electricalSize: Returns ka, the wavenumber times the radius of the smallest sphere around the plate.

            electricalSize = obj.wavenumber() * obj.mesh_.radius();
        end

        function density = density(obj)
        % density: Returns the faces per wavelength along the longer face side, the coarser direction of the mesh.

            [sizeX, sizeY] = obj.mesh_.faceSize();
            density = 2 * pi / (obj.wavenumber() * max(sizeX, sizeY));
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
                obj.requireAssembly_();
                obj.generators_ = rooftopMoM.assembler.greensMatrix.MomentAssembler().generate(obj.context_);
                rooftopMoM.resources.Budget.instance().hold(obj.generators_.memory());
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

    methods (Static)
        function problem = fromElectricalSize(electricalSize, options)
        % fromElectricalSize: Creates the problem of a rectangular plate of the electrical size ka at the wavelength
        %   1 m, k = 2 * pi, meshed with the fewest faces that reach the density along both sides.
        %
        %   Input:
        %       electricalSize: ka, with a the radius of the smallest sphere around the plate, half its diagonal.
        %
        %   Optional name-value inputs:
        %       'density': Faces per wavelength along both sides, at least (default 10).
        %       'aspect': Side along x over the side along y (default 1, a square).
        %       'quadratureOrder', 'chunkSize': As for the constructor.
        %
        %   Output:
        %       problem: rooftopMoM.Problem.

            arguments
                electricalSize (1, 1) {mustBeReal, mustBeFinite, mustBePositive}
                options.density (1, 1) {mustBeReal, mustBeFinite, mustBePositive} = 10
                options.aspect (1, 1) {mustBeReal, mustBeFinite, mustBePositive} = 1
                options.quadratureOrder {mustBeInteger, mustBePositive, mustBeOneOrTwo_} = [16, 8]
                options.chunkSize {mustBeScalarOrEmpty, mustBeReal, mustBeFinite, mustBePositive} = []
            end

            % At the wavelength 1 m the diagonal is 2 * a = ka / pi.
            aspect = double(options.aspect);
            sides = double(electricalSize) / pi * [aspect, 1] / hypot(aspect, 1);

            % Rounding must not add a face to a side that is a whole number of faces long.
            counts = max(1, ceil((1 - 4 * eps) * sides * double(options.density)));

            if all(counts == 1)
                error('rooftopMoM:problem:singleFace', ['At ka = %g and %g faces per wavelength the plate is a ' ...
                    'single face, which carries no rooftop; give a larger ''density''.'], electricalSize, ...
                    options.density);
            end

            mesh = rooftopMoM.Mesh(counts(1), counts(2), sides(1) / counts(1), sides(2) / counts(2));
            problem = rooftopMoM.Problem(mesh, 2 * pi, 'quadratureOrder', options.quadratureOrder, ...
                'chunkSize', options.chunkSize);
        end
    end

    methods (Access = private)
        function requireAssembly_(obj)
        % requireAssembly_: Raises an error if the generators and one chunk of their assembly do not fit the memory.

            generatorBytes = rooftopMoM.resources.estimateGenerators(obj.rooftops());
            chunkBytes = obj.context_.chunkSize();
            toText = @rooftopMoM.resources.utils.formatBytes;

            rooftopMoM.resources.Budget.instance().require(generatorBytes + chunkBytes, sprintf(['The assembly of %d ' ...
                'rooftops (generators %s, chunk %s)'], obj.rooftopCount(), toText(generatorBytes), toText(chunkBytes)));
        end
    end
end

function mustBeOneOrTwo_(value)
    if ~ismember(numel(value), [1, 2])
        error('rooftopMoM:problem:invalidQuadratureOrder', 'The quadrature order is one order or a pair [near, far].');
    end
end
