classdef MagneticPartAssembler < rooftopMoM.core.greensMatrix.Assembler
    properties (Access = private)
        assemblerChunkPlanner_;
    end

    methods
        function obj = MagneticPartAssembler(assemblerChunkPlanner)
            obj.assemblerChunkPlanner_ = assemblerChunkPlanner;
        end
    end

    methods
        function [Gxx, Gyy] = generateFull(obj, context)
            Gxx = obj.generateDir_(context, 'x');
            Gyy = obj.generateDir_(context, 'y');
        end

        function [Gxx, Gyy] = generateSpecified(obj, context, chunks)
            % not implemented yet
        end
    end

    methods (Access = private)
        function G = generateDir_(obj, context, direction)
            context.rooftops().countDir(direction)
        end
    end
end
