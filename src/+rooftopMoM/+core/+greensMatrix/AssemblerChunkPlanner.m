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
        function chunks = generateFull(obj, context)
        % generateFull: Returns the chunks for the chunks for electric part greens operator matrix assembly process.
            count = context.rooftops().count();

            chunks = obj.generateImpl_(context, count);
        end

        function chunks = generateFullDir(obj, context, direction)
        % generateFullDir: Returns the chunks for greens operator. Basically it utilizes the impedance matrix chunk planner because the greens matrix has the same size as the impedance matrix column.

            count = context.rooftops().countDir(direction);

            chunks = obj.generateImpl_(context, count);
        end
    end

    methods (Access = private)
        function chunks = generateImpl_(obj, context, count)
            interactionIndices = (1:count).';

            internalChunkPlanner = rooftopMoM.core.impedanceMatrix.AssemblerChunkPlanner(obj.elementCount_);

            chunks = internalChunkPlanner.generateSpecified(context, interactionIndices, count);
        end
    end
end
