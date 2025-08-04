classdef MemoryManager < handle
    properties (Access = private)
        maxMemory_ rooftopMoM.Memory;

        classHandlers_ cell;
    end

    methods
        function obj = MemoryManager(memory)
            obj.maxMemory_ = memory;
        end
    end

    methods
        function register(obj, classHandle)
            % Check if the provided class has memory method implemented.
            if ~ismethod(classHandle, 'memory')
                error('Provided class does not have `memory` method implemented.\n');
            end

            obj.classHandlers_{end + 1} = classHandle;
        end

        function memory = memoryForAssemblerGenerateMethod(context)
        end
    end
end
