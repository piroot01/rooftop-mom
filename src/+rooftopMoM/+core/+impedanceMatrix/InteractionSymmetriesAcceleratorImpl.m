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

        function [mask, interactionIndices] = generateInteractionsSingle(obj, indicesX, indicesY, gridSize, squareGrid)
            basisSize = size(indicesX, 2) + size(indicesY, 2);

            if any(gridSize == 1)
                % If any of the grid sizes is 1, then the impedance matrix is purely Toeplitz
                % SEE: Gustafsson, Mats & Tayli, Doruk & Ehrenborg, Casimir & Cismasu, Marius & Nordebo, Sven. (2016). Antenna Current Optimization using MATLAB and CVX.
                [mask, interactionIndices] = rooftopMoM.core.impedanceMatrix.InteractionSymmetriesAcceleratorImpl.generateToeplitz_(basisSize);
            else
                % Because the XY block of the impedance matrix is transpose of YX we generate it once
                [maskXY, interactionIndicesXY] = rooftopMoM.core.impedanceMatrix.InteractionSymmetriesAcceleratorImpl.generateXY_(indicesX, indicesY, gridSize);

                if squareGrid
                    [maskXX, maskYY, interactionIndicesXX, interactionIndicesYY] = rooftopMoM.core.impedanceMatrix.InteractionSymmetriesAcceleratorImpl.generateSameDirSquareGrid_(indicesX, indicesY, gridSize);
                else
                    [maskXX, maskYY, interactionIndicesXX, interactionIndicesYY] = rooftopMoM.core.impedanceMatrix.InteractionSymmetriesAcceleratorImpl.generateSameDir_(indicesX, indicesY, gridSize);
                end

                [mask, interactionIndices] = rooftopMoM.core.impedanceMatrix.InteractionSymmetriesAcceleratorImpl.generateFull_(maskXX, maskXY, maskYY, interactionIndicesXX, interactionIndicesXY, interactionIndicesYY, squareGrid);
            end
        end
    end

    methods (Access = private, Static)
        function [mask, interactionIndices] = generateToeplitz_(count)
            interactionIndices = (1:count).';
            mask = toeplitz(interactionIndices);
        end

        function [mask, interactionIndices] = generateFull_(maskXX, maskXY, maskYY, interactionIndicesXX, interactionIndicesXY, interactionIndicesYY, squareGrid)

            % Shift the indices if the maskXX and maskYY does not share any indices.
            if ~squareGrid
                maskYY = rooftopMoM.core.impedanceMatrix.utils.shiftIndices(maskYY, max(maskXX(:)));
            end

            mask = rooftopMoM.core.impedanceMatrix.utils.mergeMasks(maskXX, maskXY, maskYY);

            interactionIndices = rooftopMoM.core.impedanceMatrix.utils.mergeInteractionIndices( ...
                interactionIndicesXX, ...
                interactionIndicesXY, ...
                interactionIndicesYY, ...
                size(mask), ...
                size(maskXY), ...
                size(maskYY), ...
                size(maskXX, 1) ...
            );
        end

        function [maskXX, maskYY, interactionIndicesXX, interactionIndicesYY] = generateSameDir_(indicesX, indicesY, gridSize)
            [maskXX, interactionIndicesXX] = rooftopMoM.core.impedanceMatrix.InteractionSymmetriesAcceleratorImpl.generateSameDirSingle_(indicesX, gridSize);
            [maskYY, interactionIndicesYY] = rooftopMoM.core.impedanceMatrix.InteractionSymmetriesAcceleratorImpl.generateSameDirSingle_(indicesY, gridSize);
        end

        function [mask, interactionIndices] = generateSameDirSingle_(indices, gridSize)
        % generateSameDirSingle_: Build interaction mask and indices for same-directional basis functions.
        %
        %   Inputs:
        %       indices: 2xN array of grid linear indices for the basis function face indices.
        %       gridSize: [nX, nY] specifying the size of the full rectangular grid.
        %
        %   Outputs:
        %     mask: N x M uint32 matrix labeling unique same-directional interactions.
        %     interactionIndices: Kx2 array of patch-pair indices corresponding to nonzero mask entries.
        
            % Extract the 'positive faces' indices.
            indicesPos = indices(1, :);

            % Convert linear indices to subscripts.
            [xs, ys] = ind2sub(gridSize, indicesPos);

            % Pairwise distances in X and Y.
            dx = abs(xs - xs.');
            dy = abs(ys - ys.');

            % Combine and find unique distance‐pairs, preserving order.
            distances = [dx(:), dy(:)];
            [~, ~, indices] = unique(distances, 'rows', 'stable');

            % Reshape back into a mask.
            n = size(indices, 2);
            mask = reshape(indices, n, n);

            % Cast to uint32.
            mask = uint32(mask);

            % Interactions are just the 1…n labels.
            interactionIndices = uint32((1:n).');
        end

        function [maskXX, maskYY, interactionIndicesXX, interactionIndicesYY] = generateSameDirSquareGrid_(indicesX, indicesY, gridSize)
        % generateSameDirSquareGrid_: Build interaction mask and indices for same-directional basis functions on a square-face grid.
        %
        %   Inputs:
        %     indicesX: 2xN array of grid linear indices for horizontal faces.
        %     indicesY: 2xM array of grid linear indices for vertical faces.
        %     gridSize: [nX, nY] specifying the size of the full rectangular grid.
        %
        %   Outputs:
        %     maskXX: N x M uint32 matrix labeling unique same-directional interactions aligned with X.
        %     maskYY: N x M uint32 matrix labeling unique same-directional interactions aligned with Y.
        %     interactionIndicesXX: Kx2 array of patch-pair indices corresponding to nonzero maskXX entries.
        %     interactionIndicesYY: Kx2 array of patch-pair indices corresponding to nonzero maskYY entries.

            % Extract positive-side face indices.
            posX = indicesX(1, :);
            posY = indicesY(1, :);

            countX = size(indicesX, 2);
            countY = size(indicesY, 2);

            % Convert linear indices to subscripts.
            [xsX, ysX] = ind2sub(gridSize, posX);
            [xsY, ysY] = ind2sub(gridSize, posY);

            % Compute the pairwise distances.
            dxX = abs(xsX - xsX.');
            dyX = abs(ysX - ysX.');

            dxY = abs(xsY - xsY.');
            dyY = abs(ysY - ysY.');

            % Stack all distances.
            dX = [dxX(:), dyX(:)];
            dY = [dxY(:), dyY(:)];
            allD = [dX; dY];

            % Unique indices over combined set.
            [~, ~, allIdx] = unique(allD, 'rows', 'stable');

            % Split back
            dXCount = numel(dxX);
            idxX = allIdx(1:dXCount);
            idxY = allIdx((dXCount + 1):end);

            % Reshape to masks
            maskXX = reshape(idxX, countX, countX);
            maskYY = reshape(idxY, countY, countY);

            % Cast to uint32.
            maskXX = uint32(maskXX);
            maskYY = uint32(maskYY);

            % Unique interaction labels (first column of each mask).
            interactionIndicesXX = maskXX(:, 1);
            interactionIndicesYY_all = maskYY(:, 1);

            % Determine Y interactions not already in X.
            maskShared = ismember(interactionIndicesYY_all, interactionIndicesXX);
            interactionIndicesYY = find(~maskShared);
        end

        function [mask, interactionIndices] = generateXY_(indicesX, indicesY, gridSize)
        % generateCrossDir_: Build interaction mask and indices for cross-directional basis functions.
        %
        %   Inputs:
        %     indicesX: 2xN array of grid linear indices for horizontal faces.
        %     indicesY: 2xM array of grid linear indices for vertical faces.
        %     gridSize: [nX, nY] specifying the size of the full rectangular grid.
        %
        %   Outputs:
        %     mask: N x M uint32 matrix labeling unique cross-directional interactions.
        %     interactionIndices: Kx2 array of patch-pair indices corresponding to nonzero mask entries.

            % Convert linear indices to 2D subscripts for every cell in the grid.
            numCells = prod(gridSize);
            [allXs, allYs] = ind2sub(gridSize, 1:numCells);

            % Extract positive-side face indices for horizontal (X) and vertical (Y) faces.
            posX = indicesX(1, :);
            posY = indicesY(1, :);

            % Determine grid coordinates of each positive-side patch.
            % horz_x1, horz_y1: location of horizontal faces on the positive side.
            % horz_x2: location of the opposite (negative) side of the same face = x1 + 1.
            horz_x1 = allXs(posX).';
            horz_y1 = allYs(posX).';
            horz_x2 = horz_x1 + 1;

            % Similarly, get coordinates of vertical faces and their opposite side.
            vert_x1 = allXs(posY).';
            vert_y1 = allYs(posY).';
            vert_y2 = vert_y1 + 1;

            countX = size(indicesX, 2);
            countY = size(indicesY, 2);

            repVertX = repmat(vert_x1.', countX, 1);
            repHorzY = repmat(horz_y1, 1, countY);

            % Compute distance from horizontal face endpoints to each vertical-face coordinate.
            dist1 = abs(horz_x1 - repVertX);
            dist2 = abs(horz_x2 - repVertX);

            % Choose the further-away x-coordinate between the two horizontal endpoints.
            x_far = horz_x1 .* (dist1 > dist2) + horz_x2 .* (dist1 <= dist2);

            % Compute distances for y-axis similarly, picking the furthest y-coordinate.
            dist3 = abs(repHorzY - vert_y1.');
            dist4 = abs(repHorzY - vert_y2.');
            y_far = vert_y1.' .* (dist3 > dist4) + vert_y2.' .* (dist3 <= dist4);

            % For each face-pair, four corners of the bounding rectangle:
            % f1: intersection at (vert_x1, horz_y1),
            % f2: opposite corner at (x_far, y_far),
            % f3: other diagonal (vert_x1, y_far),
            % f4: other diagonal (x_far, horz_y1).
            f1 = sub2ind(gridSize, repVertX, repHorzY);
            f2 = sub2ind(gridSize, x_far, y_far);
            f3 = sub2ind(gridSize, repVertX, y_far);
            f4 = sub2ind(gridSize, x_far, repHorzY);

            % Stack corners into an array: size Cx4.
            corners = [f1(:), f2(:), f3(:), f4(:)];

            % Apply translational symmetry: shift all sets so their smallest index starts at 1.
            minCorner = min(corners, [], 2);
            translated = corners - (minCorner - 1) .* (minCorner ~= 1);

            % Sort pairs within each face of the rectangle for canonical ordering.
            sortedPairs = [min(translated(:,1:2), [], 2), max(translated(:,1:2), [], 2), ...
                           min(translated(:,3:4), [], 2), max(translated(:,3:4), [], 2)];

            % Identify unique interactions (rows) and produce a label mask.
            [~, ~, labelVec] = unique(sortedPairs, 'rows', 'stable');

            % Reshape label vector into a matrix of size countX x countY, convert to uint32.
            mask = uint32(reshape(labelVec, countX, countY));

            % Generate full interaction index list for nonzero mask entries.
            interactionIndices = rooftopMoM.core.impedanceMatrix.utils.generateInteractionIndicesFromMask(mask);
        end
    end
end
