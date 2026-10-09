classdef IntegrateMomentsTest < unit.Base
% IntegrateMomentsTest: Tests of rooftopMoM.assembler.greensMatrix.utils.integrateMoments.

    methods (Test)
        function testMatchesTheDefinitionOnSeparatedCells(testCase)
            % On cells three apart the Green's function is smooth, and the tensor Gauss-Legendre rule on both cells
            % integrates the definition of the moments as accurately as the reduced rule.
            [sizes, k] = deal([0.1, 0.15], 2 * pi);
            offsets = [4, 0, 3; 1, 5, 4];
            rule = rooftopMoM.assembler.greensMatrix.utils.generateReducedRule(8, [1, 1]);

            moments = rooftopMoM.assembler.greensMatrix.utils.integrateMoments(rule, offsets, sizes, k);

            [points, weights] = rooftopMoM.basisFunctions.utils.generateGaussLegendre1D(10, 0, 1);
            [xi, xiOther, eta, etaOther] = ndgrid(points, points, points, points);
            pairWeights = weights .* weights.' .* reshape(weights, 1, 1, []) .* reshape(weights, 1, 1, 1, []);

            for c = 1:size(offsets, 2)
                distance = sqrt((sizes(1) * (xi - xiOther - offsets(1, c))).^2 ...
                    + (sizes(2) * (eta - etaOther - offsets(2, c))).^2);
                values = pairWeights .* exp(-1i * k * distance) ./ distance;
                expected = [sum(values, 'all'); sum(xi .* values, 'all'); sum(xiOther .* values, 'all'); ...
                    sum(xi .* xiOther .* values, 'all'); sum(eta .* values, 'all'); sum(etaOther .* values, 'all'); ...
                    sum(eta .* etaOther .* values, 'all')];

                testCase.verifyEqual(moments(:, c), expected, 'RelTol', 1e-13);
            end
        end

        function testTouchingPairsIntegrateTheRemainder(testCase)
            [sizes, k] = deal([0.1, 0.15], 2 * pi);
            rule = rooftopMoM.assembler.greensMatrix.utils.generateReducedRule(6, [2, 1]);
            [nu, eta] = ndgrid(rule.nodesX, rule.nodesY);

            for offset = [0, 1, 0, 1; 0, 0, 1, 1]
                distance = sqrt((sizes(1) * (nu(:) - offset(1))).^2 + (sizes(2) * (eta(:) - offset(2))).^2);
                expected = rule.weights * rooftopMoM.assembler.singularIntegrals.evaluateRemainder(distance, k);

                testCase.verifyEqual(rooftopMoM.assembler.greensMatrix.utils.integrateMoments(rule, offset, sizes, k), ...
                    expected, 'RelTol', 1e-14);
            end
        end

        function testOffsetsAreIntegratedIndependently(testCase)
            rule = rooftopMoM.assembler.greensMatrix.utils.generateReducedRule(4, [1, 1]);
            offsets = [0, 1, 2, 7; 0, 1, 0, 3];
            together = rooftopMoM.assembler.greensMatrix.utils.integrateMoments(rule, offsets, [0.1, 0.1], 2 * pi);

            for c = 1:size(offsets, 2)
                alone = rooftopMoM.assembler.greensMatrix.utils.integrateMoments(rule, offsets(:, c), [0.1, 0.1], ...
                    2 * pi);
                testCase.verifyEqual(together(:, c), alone, 'RelTol', 1e-15);
            end
        end
    end
end
