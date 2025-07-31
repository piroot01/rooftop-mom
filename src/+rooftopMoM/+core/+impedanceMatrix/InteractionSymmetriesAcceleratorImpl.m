classdef InteractionSymmetriesAcceleratorImpl < handle
    methods (Access = protected)
        function [mask, interactionIndices] = generateInteractions(obj, context)
            rooftops = context.rooftops();
            mesh = rooftops.mesh();
            indices = rooftops.indices();

            gridSizes = rooftopMoM.mesh.utils.generateMeshGridSizes(mesh);

            % Determine if the mesh grid is square one.
            squareGrid = rooftopMoM.mesh.utils.isMeshFaceSquare(mesh);

            meshCount = mesh.count();
            mask = cell(1, meshCount);
            interactionIndices = cell(1, meshCount);

            for n = 1:meshCount
                [mask{n}, interactionIndices{n}] = obj.generateInteractionsSingle(indices{n}, gridSizes{n}, squareGrid(n));
            end
        end

        function [mask, interactionIndices] = generateInteractionsSingle(obj, indices, gridSizes, squareGrid)
            indices
            gridSizes
            squareGrid


            mask = [];
            interactionIndices = [];
        end
    end
end
