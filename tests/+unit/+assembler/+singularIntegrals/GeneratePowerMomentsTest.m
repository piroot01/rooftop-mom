classdef GeneratePowerMomentsTest < unit.Base
% GeneratePowerMomentsTest: Tests of rooftopMoM.assembler.singularIntegrals.generatePowerMoments.

    properties (TestParameter)
        power = {-1, 1, 3};
    end

    methods (Test)
        function testMatchesTheExactValues(testCase)
            % tests/data/exactMoments.csv holds the moments of five cells to 30 digits, by quadrature in 40-digit
            % arithmetic (tests/data/generateExactMoments.py). In double precision the corner sums of an elongated cell
            % lose digits.
            root = fileparts(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))));
            data = readmatrix(fullfile(root, 'tests', 'data', 'exactMoments.csv'), 'NumHeaderLines', 1);
            offsets = [0, 0; 1, 0; 0, 1; 1, 1];

            for row = 1:size(data, 1)
                [a, b, exponent] = deal(data(row, 1), data(row, 2), data(row, 3));
                [~, column] = ismember(data(row, 4:5), offsets, 'rows');
                moments = rooftopMoM.assembler.singularIntegrals.generatePowerMoments(a, b, exponent);
                tolerance = 1e-13;

                if max(a, b) / min(a, b) > 4
                    tolerance = 1e-12;
                end

                testCase.verifyEqual(moments(:, column), data(row, 6:12).', 'RelTol', tolerance, ...
                    sprintf('cell %g x %g, power %d, offset (%d, %d)', a, b, exponent, data(row, 4:5)));
            end
        end

        function testCoincidentCellGivesTheSelfTermsOfTheThesis(testCase)
            % The closed self terms of 1 / R in the appendix of the bachelor thesis, over two equal halves of one rooftop,
            % over different halves and over the cell; the assembler used them before the near-field extraction.
            for sizes = [0.1, 0.15; 0.3, 0.05].'
                [a, b] = deal(sizes(1), sizes(2));
                moments = rooftopMoM.assembler.singularIntegrals.generatePowerMoments(a, b, -1);

                testCase.verifyEqual(moments(1, 1), selfTermOfTheCell_(a, b), 'RelTol', 1e-14);
                testCase.verifyEqual(moments(4, 1), selfTermOfEqualHalves_(a, b), 'RelTol', 1e-14);
                testCase.verifyEqual(moments(2, 1) - moments(4, 1), selfTermOfDifferentHalves_(a, b), 'RelTol', 1e-14);
                testCase.verifyEqual(moments(7, 1), selfTermOfEqualHalves_(b, a), 'RelTol', 1e-14);
                testCase.verifyEqual(moments(5, 1) - moments(7, 1), selfTermOfDifferentHalves_(b, a), 'RelTol', 1e-14);
            end
        end

        function testTransposedCellsExchangeTheDirections(testCase, power)
            moments = rooftopMoM.assembler.singularIntegrals.generatePowerMoments(0.1, 0.25, power);
            transposed = rooftopMoM.assembler.singularIntegrals.generatePowerMoments(0.25, 0.1, power);

            % M00 comes from the corner sums of the other direction, rounded differently.
            testCase.verifyEqual(moments(5:7, :), transposed(2:4, [1, 3, 2, 4]));
            testCase.verifyEqual(moments(1, :), transposed(1, [1, 3, 2, 4]), 'RelTol', 1e-14);
        end

        function testScalesWithThePower(testCase, power)
            moments = rooftopMoM.assembler.singularIntegrals.generatePowerMoments(0.1, 0.25, power);
            scaled = rooftopMoM.assembler.singularIntegrals.generatePowerMoments(0.3, 0.75, power);

            % The corner sums of the two cells round differently.
            testCase.verifyEqual(scaled, 3^power * moments, 'RelTol', 3e-14);
        end
    end
end

function value = selfTermOfTheCell_(a, b)
    diagonal = sqrt(a^2 + b^2);
    value = 2 * (a^3 - a^2 * diagonal + b^3 - b^2 * diagonal + 3 * a * b^2 * log((a + diagonal) / b) ...
        + 3 * a^2 * b * log((b + diagonal) / a)) / (3 * a^2 * b^2);
end

function value = selfTermOfEqualHalves_(a, b)
    diagonal = sqrt(a^2 + b^2);
    value = (12 * a^5 - 12 * a^4 * diagonal + a^2 * b^2 * (30 * b - 29 * diagonal) + 2 * b^4 * (b - diagonal) ...
        + 15 * a^3 * b * (4 * b * log((a + diagonal) / b) + 3 * a * log((b + diagonal) / a))) / (90 * a^4 * b^2);
end

function value = selfTermOfDifferentHalves_(a, b)
    diagonal = sqrt(a^2 + b^2);
    value = (18 * a^5 - 2 * b^5 - 18 * a^4 * diagonal - a^2 * b^2 * diagonal + 2 * b^4 * diagonal ...
        + 30 * a^3 * b^2 * log((a + diagonal) / b) + 45 * a^4 * b * log((b + diagonal) / a)) / (90 * a^4 * b^2);
end
