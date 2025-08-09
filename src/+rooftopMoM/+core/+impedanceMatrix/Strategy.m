classdef Strategy < handle
% Strategy: High-level factory for building impedance matrices.
%
%   Notes:
%       This class provides small, well-named convenience entry points for
%       creating different parts and flavours of the rooftop MoM impedance
%       matrix (full, electric, magnetic) in either normal or accelerated form.

    methods (Static)
        function output = generateImpedanceMatrix(mesh, settings)
            % generateImpedanceMatrix: Full (magnetic + electric) matrix.
            output = rooftopMoM.core.impedanceMatrix.Strategy.generateWith( ...
                @rooftopMoM.core.impedanceMatrix.FullAssembler, mesh, settings, false ...
            );
        end

        function output = generateImpedanceMatrixMagneticPart(mesh, settings)
            % generateImpedanceMatrixMagneticPart: Only magnetic contribution.
            output = rooftopMoM.core.impedanceMatrix.Strategy.generateWith( ...
                @rooftopMoM.core.impedanceMatrix.MagneticPartAssembler, mesh, settings, false ...
            );
        end

        function output = generateImpedanceMatrixElectricPart(mesh, settings)
            % generateImpedanceMatrixElectricPart: Only electric contribution.
            output = rooftopMoM.core.impedanceMatrix.Strategy.generateWith( ...
                @rooftopMoM.core.impedanceMatrix.ElectricPartAssembler, mesh, settings, false ...
            );
        end

        function output = generateImpedanceMatrixAccelerated(mesh, settings)
            % generateImpedanceMatrixAccelerated: Full matrix using symmetry accelerator.
            output = rooftopMoM.core.impedanceMatrix.Strategy.generateWith( ...
                @rooftopMoM.core.impedanceMatrix.FullAssembler, mesh, settings, true ...
            );
        end

        function output = generateImpedanceMatrixMagneticPartAccelerated(mesh, settings)
            % generateImpedanceMatrixMagneticPartAccelerated: Magnetic part with accelerator.
            output = rooftopMoM.core.impedanceMatrix.Strategy.generateWith( ...
                @rooftopMoM.core.impedanceMatrix.MagneticPartAssembler, mesh, settings, true ...
            );
        end

        function output = generateImpedanceMatrixElectricPartAccelerated(mesh, settings)
            % generateImpedanceMatrixElectricPartAccelerated: Electric part with accelerator.
            output = rooftopMoM.core.impedanceMatrix.Strategy.generateWith( ...
                @rooftopMoM.core.impedanceMatrix.ElectricPartAssembler, mesh, settings, true ...
            );
        end
    end

    methods (Access = private, Static)
        function output = generateWith(assemblerConstructor, mesh, settings, useAccelerator)
            % generateWith: Core pipeline used by all public factory wrappers.

            [context, planner] = rooftopMoM.core.impedanceMatrix.Strategy.createAssemblerContextAndPlanner(mesh, settings);

            % Create assembler instance (expected constructor signature: assembler(planner))
            assembler = assemblerConstructor(planner);

            if useAccelerator
                % When accelerator is requested, wrap the assembler with it.
                accelerator = rooftopMoM.core.impedanceMatrix.InteractionSymmetriesAccelerator(assembler);
                output = accelerator.generateFull(context);
            else
                % Direct assembly path.
                output = assembler.generateFull(context);
            end
        end

        function [context, planner] = createAssemblerContextAndPlanner(mesh, settings)
            % createAssemblerContextAndPlanner: Prepare quadrature and rooftops.

            % Create quadrature rule from settings.
            quadrature = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(settings.order());

            % Build rooftops collector for the mesh and quadrature.
            rooftops = rooftopMoM.basisFunctions.MeshUniformRooftops(mesh, quadrature);

            % Generate the rooftop indices and intialize the quadrature points and weights.
            rooftops.generate();

            % Assemble context used by assembler.generateFull.
            context = rooftopMoM.core.impedanceMatrix.AssemblerContext( ...
                quadrature, rooftops, settings.wavenumber() ...
            );

            % Create planner controlling chunking/parallel behaviour.
            planner = rooftopMoM.core.impedanceMatrix.AssemblerChunkPlanner( ...
                settings.elementCount() ...
            );
        end
    end
end
