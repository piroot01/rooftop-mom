classdef ComputeElPart < mom.FastImpedanceMatrix.ComputeStrategy
    methods
        function result = compute(obj, strategyContext, memory_MB)
            result = obj.computeInternal(strategyContext, memory_MB);
        end
    end

    methods (Access = private)
        function result = computeInternal(obj, strategyContext, memory_MB)
            % extract constants
            k = strategyContext.constants.k;
            Z0 = strategyContext.constants.Z0;
            gridSize = strategyContext.rooftopBasisFunctions.getGridSize() + 1;
            faceCount = prod(gridSize);
            orderSquared = strategyContext.quadrature.getOrder()^2;

            [chunks, slices] = obj.getChunks(strategyContext, faceCount, memory_MB);
            chunkCount = numel(chunks);

            singularityMask = false(faceCount, 1);
            singularityMask(1) = true;

            % needs padded quadrature points
            rowQuadraturePoints = strategyContext.rooftopBasisFunctions.getQuadraturePoints(true);

            result = zeros(faceCount, 1);

            % degug
            fprintf("Computing the electric part, chunk count: %d\n", chunkCount);

            for n = 1:chunkCount
                chunk = [chunks{n}];
                slice = [slices{n}];

                chunkContext.singularityMask = singularityMask(chunk);

                sliceRange = slice(1):slice(2);

                chunkContext.chunkSize = numel(chunk);

                chunkContext.rowQuadraturePoints = rowQuadraturePoints(:, sliceRange);
                chunkContext.columnQuadraturePoints = rowQuadraturePoints(:, 1:orderSquared);

                result(chunk, :) = obj.computeChunk(strategyContext, chunkContext);

                % for debugging
                fprintf("Chunk %d computed.\n", n);
            end

            result = reshape(result, gridSize);

            % apply singularity corrections
            singularityValuesType3 = mom.ImpedanceMatrix.singularity.type3(1, 1);
            result(1, 1) = result(1, 1) + singularityValuesType3;

            % apply constants
            result = (Z0) / (1i * 4 * pi * k) * result;
        end

        function result = computeChunk(~, strategyContext, chunkContext)
            k = strategyContext.constants.k;
            orderSquared = strategyContext.quadrature.getOrder()^2;

            rowQuadraturePoints = chunkContext.rowQuadraturePoints;
            columnQuadraturePoints = chunkContext.columnQuadraturePoints;

            bigSingularityMask = repelem(chunkContext.singularityMask, orderSquared, orderSquared);

            integralKernel = mom.ImpedanceMatrix.IntegralKernel.compute(strategyContext, ...
                rowQuadraturePoints, columnQuadraturePoints);

            distanceMatrix = integralKernel.distanceMatrix;
            distanceMatrixZero = (distanceMatrix < eps);

            greensFunction = integralKernel.greensFunction;
            greensFunction(bigSingularityMask) = greensFunction(bigSingularityMask) - 1 ./ distanceMatrix(bigSingularityMask);

            integrands = greensFunction;
            integrands(distanceMatrixZero) = (-1i) * k;

            weights = repmat(strategyContext.quadrature.getWeights(), chunkContext.chunkSize, 1);
            integrals = integrands .* weights;

            integralsReshaped = reshape(integrals, orderSquared, chunkContext.chunkSize, orderSquared, 1);
            result = sum(sum(integralsReshaped, 1), 3);
        end
    end
end

