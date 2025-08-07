classdef UniformRectPrimitive < handle
    properties (Access = private)
        countX_ uint32;
        countY_ uint32;
        sizeX_ double;
        sizeY_ double;
    end

    methods
        function obj = UniformRectPrimitive(countX, countY, sizeX, sizeY)
            obj.countX_ = countX;
            obj.countY_ = countY;
            obj.sizeX_ = sizeX;
            obj.sizeY_ = sizeY;
        end
    end

    methods
        function countX = countX(obj)
            countX = obj.countX_;
        end

        function countY = countY(obj)
            countY = obj.countY_;
        end

        function sizeX = sizeX(obj)
            sizeX = obj.sizeX_;
        end

        function sizeY = sizeY(obj)
            sizeY = obj.sizeY_;
        end

        function count = count(obj)
            count = obj.countX_ * obj.countY_;
        end

        function [rectSizeX, rectSizeY] = meshSize(obj)
            rectSizeX = obj.countX_ * obj.sizeX_;
            rectSizeY = obj.countY_ * obj.sizeY_;
        end

        function [cX, cY] = gridSize(obj)
            cX = obj.countX_;
            cY = obj.countY_;
        end
    end
end
