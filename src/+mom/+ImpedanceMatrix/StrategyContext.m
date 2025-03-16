classdef StrategyContext < handle
    properties (Access = public)
        quadrature mom.Quadrature;
        basisFunctions mom.BasisFunctions;

        constants mom.ImpedanceMatrix.Constants;
    end

    methods
        function obj = StrategyContext(quadrature, basisFunctions, k)
            obj.quadrature = quadrature;
            obj.basisFunctions = basisFunctions;
            obj.constants = mom.ImpedanceMatrix.Constants(k);
        end
    end
end