classdef Constants
    properties (Access = public, Constant)
        polarities = [0, 0; 0, 1; 1, 0; 1, 1] + 1;
        polaritiesCount = 4;

        Z0 = 4 * pi * 29.9792458;
    end

    properties (Access = public)
        k double;
    end

    methods
        function obj = Constants(k)
            obj.k = k;
        end
    end
end