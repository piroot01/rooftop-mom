classdef ComputeStrategyOptimiser < handle
    properties (Access = private)
        basisFunctions_ mom.BasisFunctions;
        gridSize_ uint32;
    end

    methods
        function obj = ComputeStrategyOptimiser(basisFunctions, gridSize)
            obj.basisFunctions_ = basisFunctions;
            obj.gridSize_ = gridSize;
        end
    end

    methods
        function [mask, interactions] = getMaskAndInteractions(obj)
            [maskXX, interactionsXX] = obj.generateMaskXX(obj.basisFunctions_.getFacesX());
            [maskXY, interactionsXY] = obj.generateMaskXY(obj.basisFunctions_.getFacesX(), obj.basisFunctions_.getFacesY());

            % TODO: decide what kind of sub strategy we use
            % for now we stick to the general case
            [mask, interactions] = obj.generateFullMask(maskXX, maskXY, interactionsXX, interactionsXY);
        end
    end

    methods (Access = private)
        function [mask, interactions] = generateFullMask(~, maskXX, maskXY, interactionsXX, interactionsXY)
            interactionsXXCount = size(interactionsXX, 2);

            maskXXLastIndex = max(maskXX(:));
            maskXYLastIndex = maskXXLastIndex + max(maskXY(:));

            maskXYShifted = maskXY + maskXXLastIndex;

            mask = [maskXX, maskXYShifted; maskXYShifted.', maskXX + maskXYLastIndex];
            interactions = [interactionsXX, interactionsXY + uint32([0; interactionsXXCount]), interactionsXX + uint32([interactionsXXCount; interactionsXXCount])];
        end

        function [mask, interactions] = generateMaskXX(obj, basisFunctionsX)
            positiveFaces = basisFunctionsX(1, :);
            basisFunctionXCount = size(basisFunctionsX, 2);

            [xs, ys] = ind2sub(obj.gridSize_, positiveFaces);

            dx = abs(xs - xs.');
            dy = abs(ys - ys.');

            distances = [dx(:), dy(:)];
            [uniqueDistances, ~, indices] = unique(distances, 'rows', 'stable');
            mask = reshape(indices, basisFunctionXCount, basisFunctionXCount);
            interactions = obj.findUniqueInteractions(mask, uniqueDistances);
        end

        function [mask, interactions] = generateMaskXY(obj, basisFunctionsX, basisFunctionsY)
            [xs, ys] = ind2sub(obj.gridSize_, 1:prod(obj.gridSize_));

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

            f1 = sub2ind(obj.gridSize_, rep_x, rep_y);
            f2 = sub2ind(obj.gridSize_, x_far, y_far);
            f3 = sub2ind(obj.gridSize_, rep_x, y_far);
            f4 = sub2ind(obj.gridSize_, x_far, rep_y);

            fs = [f1(:), f2(:), f3(:), f4(:)];
            fsMin = min(fs, [], 2);
            fsTranslated = fs - (fsMin - 1) .* (fsMin ~= 1);

            fsSorted = [min(fsTranslated(:, 1:2), [], 2), max(fsTranslated(:, 1:2), [], 2), min(fsTranslated(:, 3:4), [], 2), max(fsTranslated(:, 3:4), [], 2)];
            [uniqueDistances, ~, indices] = unique(fsSorted, 'rows', 'stable');
            mask = reshape(indices, basisFunctionXCount, basisFunctionYCount);
            interactions = obj.findUniqueInteractions(mask, uniqueDistances);
        end

        function interactions = findUniqueInteractions(~, mask, uniqueDistances)
            uniqueDistancesCount = size(uniqueDistances, 1);
            interactions = zeros(2, uniqueDistancesCount, 'uint32');

            for n = 1:uniqueDistancesCount
                [r, c] = find(mask == n, 1, 'first');
                interactions(:, n) = [r; c];
            end
        end
    end
end
