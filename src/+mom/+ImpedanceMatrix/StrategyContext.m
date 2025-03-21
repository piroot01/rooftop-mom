classdef StrategyContext < handle
    properties (Access = public)
        quadrature mom.Quadrature;
        basisFunctions mom.BasisFunctions;
        rooftopBasisFunctions mom.RooftopBasisFunctions;

        constants mom.ImpedanceMatrix.Constants;
    end

    methods
        function obj = StrategyContext(quadrature, basisFunctions, rooftopBasisFunctions, k)
            obj.quadrature = quadrature;
            obj.basisFunctions = basisFunctions;
            obj.rooftopBasisFunctions = rooftopBasisFunctions;
            obj.constants = mom.ImpedanceMatrix.Constants(k);
        end
    end
end
