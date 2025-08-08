classdef AssemblerChunkPlanner < handle
    properties (Access = private)
        elementCount_ uint32;
    end

    methods
        function obj = AssemblerChunkPlanner(elementCount)
            obj.elementCount_ = elementCount;
        end
    end

    methods
        function chunks = generateDir(obj, context, direction)
        % generateFull: Returns the chunks for greens operator. Basically it utilizes the impedance matrix chunk planner because the greens matrix has the same size as the impedance matrix column.

            count = context.rooftops().countDir(direction);
            interactions = (1:count).';

            internalChunkPlanner = rooftopMoM.core.impedanceMatrix.AssemblerChunkPlanner(obj.elementCount_);

            chunks = internalChunkPlanner.generateSpecified(context, interactions, count);
        end
    end
end
