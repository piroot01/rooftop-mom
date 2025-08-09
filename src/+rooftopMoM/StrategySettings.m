classdef StrategySettings < handle
    properties (Access = private)
        option_ rooftopMoM.StrategyOptions;

        % TODO: Use memory instead
        elementCount_;
        wavenumber_;
        order_;
    end

    methods
        function obj = StrategySettings(option)
            obj.option_ = option;
        end
    end

    % Setters
    methods
        function setElementCount(obj, elementCount)
            obj.elementCount_ = elementCount;
        end

        function setWavenumber(obj, wavenumber)
            obj.wavenumber_ = wavenumber;
        end

        function setOrder(obj, order)
            obj.order_ = order;
        end
    end

    % Getters
    methods
        function option = option(obj)
            option = obj.option_;
        end

        function elementCount = elementCount(obj)
            elementCount = obj.elementCount_;
        end

        function wavenumber = wavenumber(obj)
            wavenumber = obj.wavenumber_;
        end

        function order = order(obj)
            order = obj.order_;
        end
    end
end
