classdef Strategy < handle
    properties (Access = private)
        strategySettings_ rooftopMoM.StrategySettings;
    end

    methods
        function obj = Strategy(settings)
            obj.strategySettings_ = settings;
        end
    end

    methods
        function output = generateWithMesh(obj, mesh)
            option = obj.strategySettings_.option();

            switch option
                case rooftopMoM.StrategyOptions.IMPEDANCE_MATRIX
                    output = rooftopMoM.core.impedanceMatrix.Strategy.generateImpedanceMatrix(mesh, obj.strategySettings_);
                case rooftopMoM.StrategyOptions.IMPEDANCE_MATRIX_MAGNETIC_PART
                    output = rooftopMoM.core.impedanceMatrix.Strategy.generateImpedanceMatrixMagneticPart(mesh, obj.strategySettings_);
                case rooftopMoM.StrategyOptions.IMPEDANCE_MATRIX_ELECTRIC_PART
                    output = rooftopMoM.core.impedanceMatrix.Strategy.generateImpedanceMatrixElectricPart(mesh, obj.strategySettings_);
                case rooftopMoM.StrategyOptions.IMPEDANCE_MATRIX_ACCELERATED
                    output = rooftopMoM.core.impedanceMatrix.Strategy.generateImpedanceMatrixAccelerated(mesh, obj.strategySettings_);
                case rooftopMoM.StrategyOptions.IMPEDANCE_MATRIX_MAGNETIC_PART_ACCELERATED
                    output = rooftopMoM.core.impedanceMatrix.Strategy.generateImpedanceMatrixMagneticPartAccelerated(mesh, obj.strategySettings_);
                case rooftopMoM.StrategyOptions.IMPEDANCE_MATRIX_ELECTRIC_PART_ACCELERATED
                    output = rooftopMoM.core.impedanceMatrix.Strategy.generateImpedanceMatrixElectricPartAccelerated(mesh, obj.strategySettings_);
                otherwise
                    error("Unknown strategy option '%s'.", option);
            end 
        end
    end
end
