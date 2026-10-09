classdef GenerateReducedRuleTest < unit.Base
% GenerateReducedRuleTest: Tests of rooftopMoM.assembler.greensMatrix.utils.generateReducedRule.

    properties (TestParameter)
        partCounts = {[1, 1], [3, 2]};
    end

    methods (Test)
        function testNodesLieInsideTheParts(testCase, partCounts)
            rule = rooftopMoM.assembler.greensMatrix.utils.generateReducedRule(4, partCounts);
            ends = (-partCounts(1):partCounts(1)) / partCounts(1);

            testCase.verifySize(rule.nodesX, [8 * partCounts(1), 1]);
            testCase.verifySize(rule.nodesY, [1, 8 * partCounts(2)]);
            testCase.verifySize(rule.weights, [7, 64 * prod(partCounts)]);
            testCase.verifyGreaterThan(min(abs(rule.nodesX - ends), [], 'all'), 0);
        end

        function testMatchesTheDefinitionForPolynomialKernels(testCase, partCounts)
            % A kernel polynomial in the coordinate differences, against the tensor Gauss-Legendre rule on both cells;
            % both are exact for it.
            kernel = @(alongX, alongY) (alongX - 0.3).^3 .* (alongY + 0.2).^2 + alongX .* alongY;
            rule = rooftopMoM.assembler.greensMatrix.utils.generateReducedRule(4, partCounts);
            [nu, eta] = ndgrid(rule.nodesX, rule.nodesY);

            [points, weights] = rooftopMoM.basisFunctions.utils.generateGaussLegendre1D(4, 0, 1);
            [xi, xiOther, eta1, etaOther] = ndgrid(points, points, points, points);
            pairWeights = weights .* weights.' .* reshape(weights, 1, 1, []) .* reshape(weights, 1, 1, 1, []);
            values = pairWeights .* kernel(xi - xiOther, eta1 - etaOther);
            expected = [sum(values, 'all'); sum(xi .* values, 'all'); sum(xiOther .* values, 'all'); ...
                sum(xi .* xiOther .* values, 'all'); sum(eta1 .* values, 'all'); sum(etaOther .* values, 'all'); ...
                sum(eta1 .* etaOther .* values, 'all')];

            testCase.verifyEqual(rule.weights * kernel(nu(:), eta(:)), expected, 'RelTol', 1e-13);
        end
    end
end
