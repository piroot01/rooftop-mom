classdef IntegrateMonomialTest < unit.Base
% IntegrateMonomialTest: Tests of rooftopMoM.assembler.singularIntegrals.utils.integrateMonomial.

    properties (TestParameter)
        power = {-1, 1, 3};
    end

    methods (Test)
        function testVanishesOnTheAxes(testCase, power)
            for m = 0:3
                for n = 0:1
                    values = rooftopMoM.assembler.singularIntegrals.utils.integrateMonomial(power, m, n, ...
                        [0, 0.7, 0], [0.4, 0, 0]);
                    testCase.verifyEqual(values, [0, 0, 0]);
                end
            end
        end

        function testCornerDifferencesAreTheIntegrals(testCase, power)
            % One rectangle with a corner at the singular point and one away from it.
            rectangles = [0, 0.6, 0, 0.25; 0.3, 0.8, 0.1, 0.9];

            for r = 1:size(rectangles, 1)
                [a0, a1, b0, b1] = deal(rectangles(r, 1), rectangles(r, 2), rectangles(r, 3), rectangles(r, 4));

                for m = 0:3
                    for n = 0:1
                        G = @(a, b) rooftopMoM.assembler.singularIntegrals.utils.integrateMonomial(power, m, n, a, b);
                        expected = integral2(@(u, v) u.^m .* v.^n .* (u.^2 + v.^2).^(power / 2), a0, a1, b0, b1, ...
                            'AbsTol', 0, 'RelTol', 1e-13);

                        testCase.verifyEqual(G(a1, b1) - G(a0, b1) - G(a1, b0) + G(a0, b0), expected, 'RelTol', 1e-11);
                    end
                end
            end
        end

        function testIsHomogeneous(testCase, power)
            for m = 0:3
                for n = 0:1
                    G = @(a, b) rooftopMoM.assembler.singularIntegrals.utils.integrateMonomial(power, m, n, a, b);

                    testCase.verifyEqual(G(0.9, 2.4), 3^(m + n + power + 2) * G(0.3, 0.8), 'RelTol', 1e-14);
                end
            end
        end

        function testExchangesTheArguments(testCase, power)
            G = @(m, n, a, b) rooftopMoM.assembler.singularIntegrals.utils.integrateMonomial(power, m, n, a, b);

            testCase.verifyEqual(G(0, 1, 0.3, 0.7), G(1, 0, 0.7, 0.3), 'RelTol', 1e-14);
            testCase.verifyEqual(G(1, 1, 0.3, 0.7), G(1, 1, 0.7, 0.3), 'RelTol', 1e-14);
        end
    end
end
