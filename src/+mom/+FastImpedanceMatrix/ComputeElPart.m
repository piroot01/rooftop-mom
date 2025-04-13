classdef ComputeElPart < mom.FastImpedanceMatrix.ComputeStrategy
    methods
        function result = compute(obj, strategyContext)
            result = obj.computeInternal(strategyContext);
        end
    end

    methods (Access = private)
        function result = computeInternal(~, strategyContext)
            % extract constants
            k = strategyContext.constants.k;
            Z0 = strategyContext.constants.Z0;
            gridSize = strategyContext.rooftopBasisFunctions.getGridSize() + 1;
            faceCount = prod(gridSize);
            orderSquared = strategyContext.quadrature.getOrder()^2;

            % the singular integral will always be the first one
            bigSingularityMask = false(faceCount * orderSquared, orderSquared);
            bigSingularityMask(1:orderSquared, :) = true;

            rowQuadraturePoints = strategyContext.rooftopBasisFunctions.getQuadraturePoints();
            columnQuadraturePoints = rowQuadraturePoints(:, 1:orderSquared);

            integralKernel = mom.ImpedanceMatrix.IntegralKernel.compute(strategyContext, ...
                rowQuadraturePoints, columnQuadraturePoints);

            distanceMatrix = integralKernel.distanceMatrix;
            distanceMatrixZero = (distanceMatrix < eps);

            greensFunction = integralKernel.greensFunction;
            greensFunction(bigSingularityMask) = greensFunction(bigSingularityMask) - 1 ./ distanceMatrix(bigSingularityMask);

            integrands = greensFunction;
            integrands(distanceMatrixZero) = (-1i) * k;

            weights = repmat(strategyContext.quadrature.getWeights(), faceCount, 1);
            integrals = integrands .* weights;

            integralsReshaped = reshape(integrals, orderSquared, faceCount, orderSquared, 1);
            integralsSummed = sum(sum(integralsReshaped, 1), 3);
            result = reshape(integralsSummed, gridSize);

            % apply singularity corrections
            singularityValuesType3 = mom.ImpedanceMatrix.singularity.type3(1, 1);
            result(1, 1) = result(1, 1) + singularityValuesType3;

            % apply constants
            result = (Z0) / (1i * 4 * pi * k) * result;
        end
    end
end
