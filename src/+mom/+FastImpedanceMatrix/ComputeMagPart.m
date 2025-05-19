classdef ComputeMagPart < mom.FastImpedanceMatrix.ComputeStrategy
    methods
        function [resultXDirection, resultYDirection] = compute(obj, strategyContext, memory_MB)
            resultXDirection = obj.computeXDirection(strategyContext, memory_MB);
            resultYDirection = obj.computeYDirection(strategyContext, memory_MB);
        end

        function result = computeXDirection(obj, strategyContext, memory_MB)
            % for x direction, singularity corrections are applied at
            % indices (1,1) for type1 and (2,1) for type2.
            result = obj.computeInDirectionInternal(strategyContext, 'x', [1, 1], [2, 1], memory_MB);
        end
        
        function result = computeYDirection(obj, strategyContext, memory_MB)
            % for y direction, singularity corrections are applied at
            % indices (1,1) for type1 and (1,2) for type2.
            result = obj.computeInDirectionInternal(strategyContext, 'y', [1, 1], [1, 2], memory_MB);
        end
    end

    methods (Access = private)
        function result = computeInDirectionInternal(obj, strategyContext, direction, singularityIndex1, singularityIndex2, memory_MB)
            % check if the given direction has any rooftop in the given direction
            gridSize = strategyContext.rooftopBasisFunctions.getGridSize();

            if direction == 'x' && gridSize(1) < 2
                result = nan;
                return;
            elseif direction == 'y' && gridSize(2) < 2
                result = nan;
                return;
            end

            basisFunctionCount = strategyContext.rooftopBasisFunctions.getCount();
            [chunks, slices] = obj.getChunks(strategyContext, basisFunctionCount, memory_MB);
            chunkCount = numel(chunks);

            % extract constants
            k = strategyContext.constants.k;
            Z0 = strategyContext.constants.Z0;
            p1 = strategyContext.constants.polarities(:, 1);
            p2 = strategyContext.constants.polarities(:, 2);
            orderSquared = strategyContext.quadrature.getOrder()^2;

            singularityMask = obj.getSingularityMask(strategyContext, direction, direction);

            rowQuadraturePoints = strategyContext.rooftopBasisFunctions.getQuadraturePointsInDirection(direction);
            columnQuadraturePoints = rowQuadraturePoints(:, 1:orderSquared, :);
            rowRooftopPoints = strategyContext.rooftopBasisFunctions.getRooftopPointsInDirection(direction);
            columnRooftopPoints = rowRooftopPoints(:, 1:orderSquared, :);

            result = zeros(basisFunctionCount, 1);

            % degug
            fprintf("Computing the magnetic part, chunk count: %d\n", chunkCount);

            for n = 1:chunkCount
                chunk = [chunks{n}];
                slice = [slices{n}];

                sliceRange = slice(1):slice(2);

                chunkContext.chunkSize = numel(chunk);

                chunkContext.singularityMask = singularityMask(chunk, 1, :);
                
                % size of order^2 so does not require chopping
                chunkContext.columnQuadraturePointsPolarity = columnQuadraturePoints(:, :, p1);
                chunkContext.columnRooftopPointsPolarity = columnRooftopPoints(:, :, p1);

                chunkContext.rowQuadraturePointsPolarity = rowQuadraturePoints(:, sliceRange, p2);
                chunkContext.rowRooftopPointsPolarity = rowRooftopPoints(:, sliceRange, p2);

                result(chunk, :) = obj.computeChunk(strategyContext, chunkContext);

                % for debugging
                fprintf("Chunk %d computed.\n", n);
            end

            result = reshape(result, gridSize);

            % apply singularity corrections
            singularityValuesType1 = mom.ImpedanceMatrix.singularity.type1(1, 1);
            singularityValuesType2 = mom.ImpedanceMatrix.singularity.type2(1, 1);
            result(singularityIndex1(1), singularityIndex1(2)) = result(singularityIndex1(1), singularityIndex1(2)) + 2 * singularityValuesType1;
            result(singularityIndex2(1), singularityIndex2(2)) = result(singularityIndex2(1), singularityIndex2(2)) + singularityValuesType2;

            % apply constants
            result = (1i * k * Z0) / (4 * pi) * result;
        end

        function result = computeChunk(~, strategyContext, chunkContext)
            k = strategyContext.constants.k;
            orderSquared = strategyContext.quadrature.getOrder()^2;

            integralKernel = mom.ImpedanceMatrix.IntegralKernel.compute(strategyContext, ...
                chunkContext.rowQuadraturePointsPolarity, chunkContext.columnQuadraturePointsPolarity);

            distanceMatrix = integralKernel.distanceMatrix;
            distanceMatrixZero = (distanceMatrix < eps);

            rooftopPointsDot = pagemtimes(permute(chunkContext.rowRooftopPointsPolarity, [2, 1, 3]), ...
                chunkContext.columnRooftopPointsPolarity);

            bigSingularityMask = repelem(chunkContext.singularityMask, orderSquared, orderSquared);

            greensFunction = integralKernel.greensFunction;
            greensFunction(bigSingularityMask) = greensFunction(bigSingularityMask) - 1 ./ distanceMatrix(bigSingularityMask);

            integrands = rooftopPointsDot .* greensFunction;
            integrands(distanceMatrixZero) = rooftopPointsDot(distanceMatrixZero) * (-1i) * k;

            weights = repmat(strategyContext.quadrature.getWeights(), chunkContext.chunkSize, 1);
            integrals = integrands .* weights;

            integralsReshaped = reshape(integrals, orderSquared, chunkContext.chunkSize, ...
                orderSquared, 1, strategyContext.constants.polaritiesCount);
            integralsSummed = sum(sum(integralsReshaped, 1), 3);
            results = reshape(integralsSummed, [chunkContext.chunkSize, 1, strategyContext.constants.polaritiesCount]);
            result = sum(results, 3);
        end
    end
end
