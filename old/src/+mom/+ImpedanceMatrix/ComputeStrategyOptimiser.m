classdef ComputeStrategyOptimiser < handle
    properties (Access = private)
        basisFunctions_ mom.BasisFunctions;
        rooftopUniqueUniqueSizes_ double;

        gridSizes_;
    end

    methods
        function obj = ComputeStrategyOptimiser(basisFunctions, rooftopUniqueUniqueSizes, gridSizes)
            obj.basisFunctions_ = basisFunctions;
            obj.rooftopUniqueUniqueSizes_ = rooftopUniqueUniqueSizes;

            obj.gridSizes_ = gridSizes;
        end
    end

    methods
        % for each mesh of the cluster evaluate the masks and interactions
        function [masks, interactions] = getAllMasksAndInteractions(obj)
            meshCount = obj.basisFunctions_.getIDCount();
            facesX = obj.basisFunctions_.getFacesXAsCell();
            facesY = obj.basisFunctions_.getFacesYAsCell();

            masks = cell(1, meshCount);
            interactions = cell(1, meshCount);

            for n = 1:meshCount
                sameSize = (obj.rooftopUniqueUniqueSizes_(n, 1) == obj.rooftopUniqueUniqueSizes_(n, 2));
                [masks{n}, interactions{n}] = obj.getMaskAndInteractions(facesX{n}, facesY{n}, obj.gridSizes_{n}, sameSize);
            end
        end

        % single mesh
        function [mask, interactions] = getMaskAndInteractions(obj, facesX, facesY, gridSize, sameSize)
            if any(gridSize == 1)
                [mask, interactions] = obj.generateToeplitzMask(size(facesX, 2) + size(facesY, 2));
            else
                [maskXY, interactionsXY] = obj.generateMaskXY(facesX, facesY, gridSize);

                if sameSize
                    [maskXX, maskYY, interactionsXX, interactionsYY] = obj.generateMaskXXAndYYSameSize(facesX, facesY, gridSize);
                else
                    [maskXX, maskYY, interactionsXX, interactionsYY] = obj.generateMaskXXAndYY(facesX, facesY, gridSize);
                end

                [mask, interactions] = obj.generateFullMask(maskXX, maskXY, maskYY, interactionsXX, interactionsXY, interactionsYY, sameSize);
            end
        end

        % apply the mask and generate the full mask
        function matrix = applyMask(~, matrixUniqueValues, mask, interactions)
            interactionValues = matrixUniqueValues(interactions);
            matrix = zeros(size(mask), 'like', interactionValues);
            nonZeroMask = (mask ~= 0);
            matrix(nonZeroMask) = interactionValues(mask(nonZeroMask));
            matrix(~nonZeroMask) = matrixUniqueValues(~nonZeroMask);
        end
    end

    methods (Access = private)
        function [mask, interactions] = generateToeplitzMask(~, basisFunctionCount)
            mask = toeplitz(1:basisFunctionCount);
            interactions = (1:basisFunctionCount).';
        end

        function [mask, interactions] = generateFullMask(~, maskXX, maskXY, maskYY, interactionsXX, interactionsXY, interactionsYY, sameSize)
            if ~sameSize
                maskYY = maskYY + max(maskXX(:));
            end

            maskYYLastIndex = max(maskYY(:));
            maskXYShifted = maskXY + maskYYLastIndex;

            mask = [maskXX, maskXYShifted; ...
                    maskXYShifted.', maskYY];

            function indices = transformIndices(maskSize, blockSize, blockIndices, rowStart, colStart)
                [rows, cols] = ind2sub(blockSize, blockIndices);
                rows = rows + rowStart;
                cols = cols + colStart;
                indices = uint32(sub2ind(maskSize, rows, cols));
            end

            interactionsXXCount = numel(interactionsXX);

            interactionsXYTransformed = transformIndices(size(mask), size(maskXY), interactionsXY, 0, interactionsXXCount);
            interactionsYYTransformed = transformIndices(size(mask), size(maskYY), interactionsYY, interactionsXXCount, interactionsXXCount);

            interactions = [interactionsXX; interactionsYYTransformed; interactionsXYTransformed];
        end

        function [maskXX, maskYY, interactionsXX, interactionsYY] = generateMaskXXAndYY(~, basisFunctionsX, basisFunctionsY, gridSize)
            positiveFacesX = basisFunctionsX(1, :);
            positiveFacesY = basisFunctionsY(1, :);
            basisFunctionXCount = size(basisFunctionsX, 2);
            basisFunctionYCount = size(basisFunctionsY, 2);

            [xsX, ysX] = ind2sub(gridSize, positiveFacesX);
            [xsY, ysY] = ind2sub(gridSize, positiveFacesY);

            dxX = abs(xsX - xsX.');
            dyX = abs(ysX - ysX.');

            dxY = abs(xsY - xsY.');
            dyY = abs(ysY - ysY.');

            distancesX = [dxX(:), dyX(:)];
            distancesY = [dxY(:), dyY(:)];

            [~, ~, indicesX] = unique(distancesX, 'rows', 'stable');
            [~, ~, indicesY] = unique(distancesY, 'rows', 'stable');

            maskXX = reshape(indicesX, basisFunctionXCount, basisFunctionXCount);
            maskYY = reshape(indicesY, basisFunctionYCount, basisFunctionYCount);

            interactionsXX = (1:size(maskXX, 1)).';
            interactionsYY = (1:size(maskYY, 1)).';
        end

        function [maskXX, maskYY, interactionsXX, interactionsYY] = generateMaskXXAndYYSameSize(~, basisFunctionsX, basisFunctionsY, gridSize)
            positiveFacesX = basisFunctionsX(1, :);
            positiveFacesY = basisFunctionsY(1, :);
            basisFunctionXCount = size(basisFunctionsX, 2);
            basisFunctionYCount = size(basisFunctionsY, 2);

            [xsX, ysX] = ind2sub(gridSize, positiveFacesX);
            [ysY, xsY] = ind2sub(gridSize, positiveFacesY);

            dxX = abs(xsX - xsX.');
            dyX = abs(ysX - ysX.');

            dxY = abs(xsY - xsY.');
            dyY = abs(ysY - ysY.');

            distancesX = [dxX(:), dyX(:)];
            distancesY = [dxY(:), dyY(:)];
            distances = [distancesX; distancesY];

            [~, ~, indices] = unique(distances, 'rows', 'stable');

            indicesXCount = numel(dxX);
            indicesX = indices(1:indicesXCount);
            indicesY = indices((indicesXCount + 1):end);

            maskXX = reshape(indicesX, basisFunctionXCount, basisFunctionXCount);
            maskYY = reshape(indicesY, basisFunctionYCount, basisFunctionYCount);

            maskXXUniqueInteractions = maskXX(:, 1);
            maskYYUniqueInteractions = maskYY(:, 1);

            indicesX = ismember(maskYYUniqueInteractions, maskXXUniqueInteractions);

            interactionsXX = maskXXUniqueInteractions;
            interactionsYY = find(~indicesX);
        end

        function [mask, interactions] = generateMaskXY(obj, basisFunctionsX, basisFunctionsY, gridSize)
            [xs, ys] = ind2sub(gridSize, 1:prod(gridSize));

            horz_x1 = xs(basisFunctionsX(1, :)).';
            horz_x2 = xs(basisFunctionsX(2, :)).';
            horz_y1 = ys(basisFunctionsX(1, :)).';

            vert_x1 = xs(basisFunctionsY(1, :)).';
            vert_y1 = ys(basisFunctionsY(1, :)).';
            vert_y2 = ys(basisFunctionsY(2, :)).';

            [basisFunctionXCount, basisFunctionYCount] = deal(size(basisFunctionsX, 2), size(basisFunctionsY, 2));
            rep_x = repmat(vert_x1.', basisFunctionXCount, 1);
            rep_y = repmat(horz_y1, 1, basisFunctionYCount);

            d1 = abs(horz_x1 - rep_x);
            d2 = abs(horz_x2 - rep_x);
            x_far = horz_x1 .* (d1 > d2) + horz_x2 .* (d1 <= d2);

            d3 = abs(rep_y - vert_y1.');
            d4 = abs(rep_y - vert_y2.');
            y_far = vert_y1.' .* (d3 > d4) + vert_y2.' .* (d3 <= d4);

            f1 = sub2ind(gridSize, rep_x, rep_y);
            f2 = sub2ind(gridSize, x_far, y_far);
            f3 = sub2ind(gridSize, rep_x, y_far);
            f4 = sub2ind(gridSize, x_far, rep_y);

            fs = [f1(:), f2(:), f3(:), f4(:)];
            fsMin = min(fs, [], 2);
            fsTranslated = fs - (fsMin - 1) .* (fsMin ~= 1);

            fsSorted = [min(fsTranslated(:, 1:2), [], 2), max(fsTranslated(:, 1:2), [], 2), min(fsTranslated(:, 3:4), [], 2), max(fsTranslated(:, 3:4), [], 2)];
            [~, ~, indices] = unique(fsSorted, 'rows', 'stable');
            mask = reshape(indices, basisFunctionXCount, basisFunctionYCount);
            interactions = obj.findUniqueInteractions(mask);
        end

        function interactions = findUniqueInteractions(~, mask)
            interactions = uint32(accumarray(mask(:), (1:numel(mask)).', [], @min));
        end
    end
end
