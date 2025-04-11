classdef ComputeMagPart < mom.FastImpedanceMatrix.ComputeStrategy
    methods
        function result = computeXDirection(obj, strategyContext)
            k = strategyContext.constants.k;
            Z0 = strategyContext.constants.Z0;

            p1 = strategyContext.constants.polarities(:, 1);
            p2 = strategyContext.constants.polarities(:, 2);

            basisFunctionCount = strategyContext.rooftopBasisFunctions.getCount();
            orderSquared = strategyContext.quadrature.getOrder()^2;

            singularityMask = obj.getSingularityMask(strategyContext, 'x', 'x');

            rowQuadraturePoints = strategyContext.rooftopBasisFunctions.getQuadraturePointsInDirection('x');
            columnQuadraturePoints = rowQuadraturePoints(:, 1:orderSquared, :);

            rowRooftopPoints = strategyContext.rooftopBasisFunctions.getRooftopPointsInDirection('x');
            columnRooftopPoints = rowRooftopPoints(:, 1:orderSquared, :);

            columnQuadraturePointsPolarity = columnQuadraturePoints(:, :, p1);
            rowQuadraturePointsPolarity = rowQuadraturePoints(:, :, p2);

            columnRooftopPointsPolarity = columnRooftopPoints(:, :, p1);
            rowRooftopPointsPolarity = rowRooftopPoints(:, :, p2);

            integralKernel = mom.ImpedanceMatrix.IntegralKernel.compute(strategyContext, rowQuadraturePointsPolarity, columnQuadraturePointsPolarity);

            distanceMatrix = integralKernel.distanceMatrix;
            distanceMatrixZero = (distanceMatrix < eps);

            rooftopPointsDot = pagemtimes(permute(rowRooftopPointsPolarity, [2, 1, 3]), columnRooftopPointsPolarity);

            bigSingularityMask = repelem(singularityMask, orderSquared, orderSquared);

            greensFunction = integralKernel.greensFunction;
            greensFunction(bigSingularityMask) = greensFunction(bigSingularityMask) - 1 ./ distanceMatrix(bigSingularityMask);

            integrands = rooftopPointsDot .* greensFunction;
            integrands(distanceMatrixZero) = rooftopPointsDot(distanceMatrixZero) * (-1i) * k;

            weights = repmat(strategyContext.quadrature.getWeights(), basisFunctionCount, 1);

            integrals = integrands .* weights;
            integralsReshaped = reshape(integrals, orderSquared, basisFunctionCount, orderSquared, 1, strategyContext.constants.polaritiesCount);
            integralsSummed = sum(sum(integralsReshaped, 1), 3);
            results = reshape(integralsSummed, [strategyContext.rooftopBasisFunctions.getGridSize(), strategyContext.constants.polaritiesCount]);

            result = sum(results, 3);

            % apply singularities
            singularityValuesType1 = mom.ImpedanceMatrix.singularity.type1(1, 1);
            singularityValuesType2 = mom.ImpedanceMatrix.singularity.type2(1, 1);

            result(1, 1) = result(1,1) + 2 * singularityValuesType1;
            result(2, 1) = result(2, 1) + singularityValuesType2;

            % apply constants
            result = (1i * k * Z0) / (4 * pi) * result;
        end
    end
end
