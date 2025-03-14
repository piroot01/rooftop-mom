classdef StrategyContext < handle
    properties (Access = public)
        quadrature mom.Quadrature;
        basisFunctions mom.BasisFunctions;

        constants = mom.ImpedanceMatrix.Constants();
    end

    methods
        function obj = StrategyContext(quadrature, basisFunctions)
            obj.quadrature = quadrature;
            obj.basisFunctions = basisFunctions;
        end
    end
end