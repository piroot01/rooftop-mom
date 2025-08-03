classdef MemoryManager < handle
    properties (Access = private)
        maxMemoryBytes_;

        classHandler_ cell;
    end

    methods
        function obj = MemoryManager(maxMemory, unit)
            obj.maxMemoryBytes_ = rooftopMoM.core.utils.convertToBytes(maxMemory, unit);
        end
    end

    methods
        function register(obj, classHandle)
            % Check if the provided class has memory method implemented.
            if ismethod(classHandle, 'memory')
                error('Provided class does not have `memory` method implemented.\n');
            end

            obj.classHandler_{end + 1} = classHandle;
        end
    end
end
