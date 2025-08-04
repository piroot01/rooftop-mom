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
            % Basically compute the remaining memory for the matrix computation.
            % Then reduce the memory by the numel(k) * size(Z).
            % Subsequently reduce the memory by extra large arrays (potentially the masks, extracted arrays, etc.)
            % From the remaining memory compute the largest distance matrix that is possible to fit.
            % Keep in mind that the distace matrix is calculated by 4 copies at once.
        end
    end
end
