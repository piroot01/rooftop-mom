classdef Strategy < handle
    methods (Static)
        function output = generateImpedanceMatrix(mesh, settings)
            output = rooftopMoM.core.impedanceMatrix.Strategy.runAssembler( ...
                @rooftopMoM.core.impedanceMatrix.FullAssembler, ...
                mesh, settings ...
            );
        end

        function output = generateImpedanceMatrixMagneticPart(mesh, settings)
            output = rooftopMoM.core.impedanceMatrix.Strategy.runAssembler( ...
                @rooftopMoM.core.impedanceMatrix.MagneticPartAssembler, ...
                mesh, settings ...
            );
        end

        function output = generateImpedanceMatrixElectricPart(mesh, settings)
            output = rooftopMoM.core.impedanceMatrix.Strategy.runAssembler( ...
                @rooftopMoM.core.impedanceMatrix.ElectricPartAssembler, ...
                mesh, settings ...
            );
        end
    end

    methods (Access = private, Static)
        function output = runAssembler(assemblerCtor, mesh, settings)
            [context, planner] = rooftopMoM.core.impedanceMatrix.Strategy.generateAssemblerContextAndPlanner(mesh, settings);
            assembler = assemblerCtor(planner);
            output = assembler.generateFull(context);
        end

        function [context, planner] = generateAssemblerContextAndPlanner(mesh, settings)
            quadrature = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(settings.order());
            rooftops = rooftopMoM.basisFunctions.MeshUniformRooftops(mesh, quadrature);
            rooftops.generate();

            context = rooftopMoM.core.impedanceMatrix.AssemblerContext( ...
                quadrature, rooftops, settings.wavenumber() ...
            );

            planner = rooftopMoM.core.impedanceMatrix.AssemblerChunkPlanner( ...
                settings.elementCount() ...
            );
        end
    end
end

