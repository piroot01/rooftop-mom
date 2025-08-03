classdef DataTypeMemory
    enumeration
        Double
        Single
        ComplexDouble
        ComplexSingle
        Int32
        UInt32
        Int64
        UInt64
        Int16
        UInt16
        Int8
        UInt8
        Logical
    end

    methods (Static)
        function memory = memory(dataType)
            switch dataType
                case rooftopMoM.DataTypeMemory.Double
                    memory = rooftopMoM.Memory(8, 'B');
                case rooftopMoM.DataTypeMemory.Single
                    memory = rooftopMoM.Memory(4, 'B');
                case rooftopMoM.DataTypeMemory.ComplexDouble
                    memory = rooftopMoM.Memory(16, 'B');
                case rooftopMoM.DataTypeMemory.ComplexSingle
                    memory = rooftopMoM.Memory(8, 'B');
                case rooftopMoM.DataTypeMemory.Int64
                    memory = rooftopMoM.Memory(8, 'B');
                case rooftopMoM.DataTypeMemory.UInt64
                    memory = rooftopMoM.Memory(8, 'B');
                case rooftopMoM.DataTypeMemory.Int32
                    memory = rooftopMoM.Memory(4, 'B');
                case rooftopMoM.DataTypeMemory.UInt32
                    memory = rooftopMoM.Memory(4, 'B');
                case rooftopMoM.DataTypeMemory.Int16
                    memory = rooftopMoM.Memory(2, 'B');
                case rooftopMoM.DataTypeMemory.UInt16
                    memory = rooftopMoM.Memory(2, 'B');
                case rooftopMoM.DataTypeMemory.Int8
                    memory = rooftopMoM.Memory(1, 'B');
                case rooftopMoM.DataTypeMemory.UInt8
                    memory = rooftopMoM.Memory(1, 'B');
                case rooftopMoM.DataTypeMemory.Logical
                    memory = rooftopMoM.Memory(1, 'B');
                otherwise
                    error('Unknown data type.');
            end
        end
    end
end
