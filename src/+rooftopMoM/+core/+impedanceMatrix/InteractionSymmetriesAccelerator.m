classdef InteractionSymmetriesAccelerator < rooftopMoM.core.impedanceMatrix.InteractionSymmetriesAcceleratorImpl
    properties (Access = private)
        strategy_;
    end

    methods
        function obj = InteractionSymmetriesAccelerator(strategy)
            obj.strategy_ = strategy;
        end
    end

    methods
        function Z = generateFull(obj, context)
            % TODO: Add support of multiple meshes.
            [mask, interactionIndices] = obj.generateInteractions(context);

            Z = [];
        end
    end
end
