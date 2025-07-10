classdef UniformRectPrimitive < handle
    properties (Access = private)
        countX_ uint32;
        countY_ uint32;
        sizeX_ double;
        sizeY_ double;

        genom_ logical;
    end

    methods
        function obj = UniformRectPrimitive(countX, countY, sizeX, sizeY, genom)
            obj.countX_ = countX;
            obj.countY_ = countY;
            obj.sizeX_ = sizeX;
            obj.sizeY_ = sizeY;

            if ~isempty(genom) && numel(genom) == count(obj)
                obj.genom_ = logical(genom);
            else
                obj.genom_ = true(1, obj.countX_ * obj.countY_);
            end
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

        function count = countGenome(obj)
            count = sum(obj.genom_);
        end

        function indices = indices(obj)
            indices = 1:count(obj);
        end

        function indicesGenome = indicesGenome(obj)
            indices = obj.indices();
            indicesGenome = indices(obj.genom_);
        end

        function [rectSizeX, rectSizeY] = meshSize(obj)
            rectSizeX = obj.countX_ * obj.sizeX_;
            rectSizeY = obj.countY_ * obj.sizeY_;
        end
    end
end
