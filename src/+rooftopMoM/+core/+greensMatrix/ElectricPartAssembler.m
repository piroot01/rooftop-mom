classdef ElectricPartAssembler < rooftopMoM.core.greensMatrix.Assembler
    properties (Access = private)
        assemblerChunkPlanner_;
    end

    methods
        function obj = ElectricPartAssembler(assemblerChunkPlanner)
            obj.assemblerChunkPlanner_ = assemblerChunkPlanner;
        end
    end

    methods
        function G = generateFull(obj, context)
            % Produce all chunks for full matrix assembly
            chunks = obj.assemblerChunkPlanner_.generateFull(context);

            [Gxx, Gyy] = obj.generateSpecified(context, chunksX, chunksY);
        end

        function [Gxx, Gyy] = generateSpecified(obj, context, chunksX, chunksY)
            Gxx = obj.generateDirSpecified(context, chunksX, 'x');
            Gyy = obj.generateDirSpecified(context, chunksY, 'y');
        end
    end
end
