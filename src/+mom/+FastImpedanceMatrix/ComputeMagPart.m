classdef ComputeMagPart < mom.FastImpedanceMatrix.ComputeStrategy
    methods
        function result = computeXDirection(obj, strategyContext)
            k = strategyContext.constants.k;
            Z0 = strategyContext.constants.Z0;
            p1 = strategyContext.constants.polarities(:, 1);
            p2 = strategyContext.constants.polarities(:, 2);
            orderSquared = strategyContext.quadrature.getOrder()^2;
            singularityMask = obj.getSingularityMask(strategyContext, 'x', 'x', 1);

            result = zeros()
            rowQuadraturePoints = strategyContext.rooftopBasisFunctions.getQuadraturePointsInDirection('x');
            columnQuadraturePoints = rowQuadraturePoints(:, 1:orderSquared, :);

            rowRooftopPoints = strategyContext.rooftopBasisFunctions.getRooftopPointsInDirection('x');
            columnRooftopPoints = rowRooftopPoints(:, 1:orderSquared, :);

            rowQuadraturePointsPolarity = rowQuadraturePoints(:, :, p1);
            columnQuadraturePointsPolarity = columnQuadraturePoints(:, :, p2);

            rowRooftopPointsPolarity = rowRooftopPoints(:, :, p1);
            columnRooftopPointsPolarity = columnRooftopPoints(:, :, p2);

            integralKernel = mom.ImpedanceMatrix.IntegralKernel.compute(strategyContext, rowQuadraturePointsPolarity, columnQuadraturePointsPolarity)

            distanceMatrix = integralKernel.distanceMatrix;
            distanceMatrixZero = (distanceMatrix < eps);

            rooftopPointsDot = pagemtimes(permute(rowRooftopPointsPolarity, [2, 1, 3]), columnRooftopPointsPolarity);
            rooftopPointsDotZero = (rooftopPointsDot < eps);


        end
    end
end
