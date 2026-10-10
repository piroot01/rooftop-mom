classdef Problem < handle
% Problem: The electric field integral equation of a Mesh at one wavenumber, discretized with rooftops. The
%   generators of its Green's matrix are assembled once, when an Operator or a Preconditioner needs them first.
%       problem = rooftopMoM.Problem(mesh, 2 * pi, 'quadratureOrder', [16, 8]);

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
                wavenumber (1, 1) {mustBeReal, mustBeFinite, mustBePositive}
                options.quadratureOrder {mustBeInteger, mustBePositive, mustBeOneOrTwo_} = [16, 8]
                options.chunkSize {mustBeScalarOrEmpty, mustBeReal, mustBeFinite, mustBePositive} = []
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
                'rad/m, quadrature orders %d near and %d far.'], countX, countY, sizeX, sizeY, rooftops.count(), ...
                wavenumber, orders(1), orders(2));
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
