classdef UniformRectPrimitive
% UniformRectPrimitive: Uniform grid of countX x countY rectangular faces of the size sizeX x sizeY.
%   The grid lies in the plane z = 0 with one corner at the origin.

    properties (Access = private)
        countX_ double;
        countY_ double;
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
        % count: Returns the number of faces.

            count = obj.countX_ * obj.countY_;
        end
    end
end
