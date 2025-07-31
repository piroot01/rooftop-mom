classdef InteractionSymmetriesAcceleratorImpl < handle
    methods (Access = protected)
        function [mask, interactionIndices] = generateInteractions(obj, context)
            rooftops = context.rooftops();
            mesh = rooftops.mesh();
            indicesX = rooftops.indicesDir('x');
            indicesY = rooftops.indicesDir('y');

            gridSizes = rooftopMoM.mesh.utils.generateMeshGridSizes(mesh);

            % Determine if the mesh grid is square one.
            squareGrid = rooftopMoM.mesh.utils.isMeshFaceSquare(mesh);

            meshCount = mesh.count();
            mask = cell(1, meshCount);
            interactionIndices = cell(1, meshCount);

            for n = 1:meshCount
                [mask{n}, interactionIndices{n}] = obj.generateInteractionsSingle(indicesX{n}, indicesY{n}, gridSizes{n}, squareGrid(n));
            end
        end

        function [mask, interactionIndices] = generateInteractionsSingle(obj, indicesX, indicesY, gridSizes, squareGrid)
            indicesX
            indicesY
            gridSizes
            squareGrid


            mask = [];
            interactionIndices = [];
        end
    end
end
