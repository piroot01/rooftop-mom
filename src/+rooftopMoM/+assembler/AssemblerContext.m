classdef AssemblerContext
% AssemblerContext: Inputs of the Green's matrix assembly: the rooftop basis on its grid, the wavenumber, the memory
%   budget of one chunk of the assembly loops and the quadrature orders of the near and the far cell pairs.

    properties (Access = private)
        rooftops_ rooftopMoM.basisFunctions.UniformRooftops;
        wavenumber_ double;
        memoryBytes_ double;
        quadratureOrders_ double;
    end

    methods
        function obj = AssemblerContext(rooftops, wavenumber, memoryBytes, quadratureOrders)
            obj.rooftops_ = rooftops;
            obj.wavenumber_ = wavenumber;
            obj.memoryBytes_ = memoryBytes;
            obj.quadratureOrders_ = quadratureOrders;
        end
    end

    methods
        function rooftops = rooftops(obj)
            rooftops = obj.rooftops_;
        end

        function wavenumber = wavenumber(obj)
            wavenumber = obj.wavenumber_;
        end

        function memoryBytes = memoryBytes(obj)
            memoryBytes = obj.memoryBytes_;
        end

        function quadratureOrders = quadratureOrders(obj)
        % quadratureOrders: Returns the 1x2 Gauss-Legendre orders [near, far] of the cell pairs (MomentAssembler).

            quadratureOrders = obj.quadratureOrders_;
        end

        function count = chunkLength(obj, bytesPerItem)
        % chunkLength: Returns how many items of a chunked loop fit into the memory budget, at least one.
        %
        %   Input:
        %       bytesPerItem: Bytes of the temporaries of one item.
        %
        %   Output:
        %       count: Positive integer.

            count = max(1, floor(obj.memoryBytes_ / bytesPerItem));
        end
    end
end
