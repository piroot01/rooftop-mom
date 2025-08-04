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
        function registerContext(obj, context)
            obj.registerSingle(context.rooftops().mesh());
            obj.registerSingle(context.rooftops());
            obj.registerSingle(context.quadrature());
        end

        function registerSingle(obj, classHandle)
            % Check if the provided class has memory method implemented.
            if ~ismethod(classHandle, 'memory')
                error('Provided class does not have `memory` method implemented.\n');
            end

            obj.classHandlers_{end + 1} = classHandle;
        end

        function count = elementCountMagneticPart(obj, context)
            remainingMemory = obj.remainingMemoryMagneticPart_(context);

            count = rooftopMoM.core.impedanceMatrix.utils.elementCountMagneticPart(context, remainingMemory);
        end
    end

    methods (Access = private)
        function memory = remainingMemoryMagneticPart_(obj, context)
            usedMemory = rooftopMoM.Memory(0, 'B');

            for n = 1:numel(obj.classHandlers_)
                classMemory = obj.classHandlers_{n}.memory();
                usedMemory = usedMemory + classMemory;
            end

            rooftops = context.rooftops();
            count = rooftops.count();

            % Memory required for the global variables of the impedance matrix assembly procedure.
            rooftopPointsMemory = rooftopMoM.core.impedanceMatrix.utils.memoryForRooftopPoints(rooftops);
            rooftopValuesMemory = rooftopMoM.core.impedanceMatrix.utils.memoryForRooftopValues(rooftops);
            idMaskMemory = count^2 * rooftopMoM.DataTypeMemory.memory(rooftopMoM.DataTypeMemory.UInt32);
            jacobiansMemory = count * rooftopMoM.DataTypeMemory.memory(rooftopMoM.DataTypeMemory.Double);
            sizeMaskMemory = count^2 * rooftopMoM.DataTypeMemory.memory(rooftopMoM.DataTypeMemory.UInt32);
            indexMaskMemory = count^2 * 4 * rooftopMoM.DataTypeMemory.memory(rooftopMoM.DataTypeMemory.Logical);

            memory = obj.maxMemory_ - usedMemory ...
                - rooftopPointsMemory ...
                - rooftopValuesMemory ...
                - idMaskMemory ...
                - jacobiansMemory ...
                - sizeMaskMemory ...
                - indexMaskMemory;
        end
    end
end
