classdef Vec3
    properties
        x double = 0;
        y double = 0;
        z double = 0;
    end

    methods
        function obj = Vec3(x, y, z)
            obj.x = x;
            obj.y = y;
            obj.z = z;
        end
    end

    methods
        function v = getAsColVec(obj)
            v = [obj.x; obj.y; obj.z];
        end

        function disp(obj)
            fprintf('[%g, %g, %g]', obj.x, obj.y, obj.z);
        end
    end
end