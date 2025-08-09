classdef Strategy < handle
    properties (Access = private)
        strategySettings_ rooftopMoM.StrategySettings;
    end

    methods
        function obj = Strategy(strategySettings)
            obj.strategySettings_ = strategySettings;
        end
    end

    methods
        function strategyOutput = generate(obj)
        end
    end
end
