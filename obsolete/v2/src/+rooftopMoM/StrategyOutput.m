classdef StrategyOutput < handle
    properties (Access = protected)
        value_;
    end

    methods
        function obj = StrategyOutput(value)
            obj.value_ = value;
        end
    end

    methods
        function setValue(obj, value)
            obj.value_ = value;
        end
    end

    methods
        function value = value(obj)
            value = obj.value_;
        end
    end
end
