classdef FullAssembler < rooftopMoM.core.greensMatrix.Assembler
    properties (Access = private)
        assemblerChunkPlanner_;
    end

    methods
        function obj = FullAssembler(assemblerChunkPlanner)
            obj.assemblerChunkPlanner_ = assemblerChunkPlanner;
        end
    end

    methods
        function output = generateFull(obj, context)
            primitive = context.rooftops().primitive();

            magneticPartAssembler = rooftopMoM.core.greensMatrix.MagneticPartAssembler(obj.assemblerChunkPlanner_);
            electricPartAssembler = rooftopMoM.core.greensMatrix.ElectricPartAssembler(obj.assemblerChunkPlanner_);

            % TODO: When the size of the grid in at leat one direction is 1 then the magnetic part assembly process crashes.
            [GmXX, GmYY] = magneticPartAssembler.generateFull(context)
            Ge = electricPartAssembler.generateFull(context)

            [M, N] = primitive.gridSize();

            % Perform the toeplitz kernel embedding into the kernel of the circulant tensor.
            GmXX_circ = rooftopMoM.core.greensMatrix.utils.embedToeplitz4ToCirculant(GmXX, M, N)
            GmYY_circ = rooftopMoM.core.greensMatrix.utils.embedToeplitz4ToCirculant(GmYY, M, N)
            Ge_circ = rooftopMoM.core.greensMatrix.utils.embedToeplitz4ToCirculant(Ge, M, N)

            output = [];
        end

        function generateSpecified(~)
            % NOTE: not implemented
        end
    end
end
