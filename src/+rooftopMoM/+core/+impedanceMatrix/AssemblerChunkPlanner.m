classdef AssemblerChunkPlanner < handle
    properties (Access = private)
        doubleCount_ uint32;
    end

    methods
        function obj = AssemblerChunkPlanner(doubleCount)
            obj.doubleCount_ = doubleCount;
        end
    end

    methods
        function chunks = generateFull(rooftops, includeDiagonalBlocks)
        
        end
    
        function chunks = generateSpecified(rooftops, interactions)
        end
    end
end
