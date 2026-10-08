classdef IntegrateFacePairsTest < unit.Base
% IntegrateFacePairsTest: Tests of rooftopMoM.assembler.greensMatrix.utils.integrateFacePairs.

    methods (Test)
        function testMatchesAnExplicitSum(testCase)
            context = testCase.createContext(4, 3, 0.1, 0.15, 3, 2 * pi);
            weights = context.rooftops().quadrature().weights4D();
            shift = [1; -1];

            integrals = rooftopMoM.assembler.greensMatrix.utils.integrateFacePairs(context, [4, 3], shift, weights, 0);

            for i = 0:3
                for j = 0:2
                    if all([i; j] + shift == 0)
                        continue;
                    end

                    expected = testCase.explicitSum_(context, [i; j] + shift, weights);
                    testCase.verifyEqual(integrals(i + 1, j + 1), expected, 'RelTol', 1e-13);
                end
            end
        end

        function testCoincidentPairGetsTheStaticIntegral(testCase)
            context = testCase.createContext(3, 3, 0.1, 0.15, 3, 2 * pi);
            weights = context.rooftops().quadrature().weights4D();

            withStatic = rooftopMoM.assembler.greensMatrix.utils.integrateFacePairs(context, [3, 3], [0; 0], weights, 5);
            withoutStatic = rooftopMoM.assembler.greensMatrix.utils.integrateFacePairs(context, [3, 3], [0; 0], weights, 0);

            testCase.verifyEqual(withStatic(1, 1) - withoutStatic(1, 1), 5, 'AbsTol', 1e-14);
            testCase.verifyEqual(withStatic(2:end), withoutStatic(2:end));
        end

        function testCoincidentPairOutsideTheGridIsSkipped(testCase)
            context = testCase.createContext(3, 3, 0.1, 0.15, 3, 2 * pi);
            weights = context.rooftops().quadrature().weights4D();

            withStatic = rooftopMoM.assembler.greensMatrix.utils.integrateFacePairs(context, [3, 3], [1; 0], weights, 5);
            withoutStatic = rooftopMoM.assembler.greensMatrix.utils.integrateFacePairs(context, [3, 3], [1; 0], weights, 0);

            testCase.verifyEqual(withStatic, withoutStatic);
        end

        function testChunkingDoesNotChangeTheResult(testCase)
            weights = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(3).weights4D();
            oneChunk = testCase.createContext(5, 4, 0.1, 0.15, 3, 2 * pi, 2^30);
            oneEntryPerChunk = testCase.createContext(5, 4, 0.1, 0.15, 3, 2 * pi, 1);

            testCase.verifyEqual( ...
                rooftopMoM.assembler.greensMatrix.utils.integrateFacePairs(oneEntryPerChunk, [5, 4], [0; 0], weights, 1), ...
                rooftopMoM.assembler.greensMatrix.utils.integrateFacePairs(oneChunk, [5, 4], [0; 0], weights, 1), ...
                'RelTol', 1e-14);
        end
    end

    methods (Access = private)
        function value = explicitSum_(~, context, offset, weights)
            primitive = context.rooftops().primitive();
            [pointsU, pointsV] = context.rooftops().quadrature().points();
            value = 0;

            for b = 1:numel(pointsU)
                for c = 1:numel(pointsU)
                    distance = hypot(primitive.sizeX() * (offset(1) + pointsU(b) - pointsU(c)), ...
                        primitive.sizeY() * (offset(2) + pointsV(b) - pointsV(c)));
                    value = value + weights(b, c) * exp(-1i * context.wavenumber() * distance) / distance;
                end
            end
        end
    end
end
