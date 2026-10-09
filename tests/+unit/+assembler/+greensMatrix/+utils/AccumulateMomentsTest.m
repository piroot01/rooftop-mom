classdef AccumulateMomentsTest < unit.Base
% AccumulateMomentsTest: Tests of rooftopMoM.assembler.greensMatrix.utils.accumulateMoments.

    methods (Test)
        function testBuildsTheGeneratorsFromTheMoments(testCase)
            % Random moments of all offsets of a 5 x 4 grid, added in two chunks in random order, against the
            % combination of the near-field note (Generators from moments) on the whole arrays.
            counts = [5, 4];
            rng(1);
            moments = complex(randn(7, prod(counts)), randn(7, prod(counts)));
            [m, n] = ind2sub(counts, 1:prod(counts));
            offsets = [m - 1; n - 1];
            parts = struct('electric', zeros(counts, 'like', 1i), 'magneticXX', zeros(counts - [1, 0], 'like', 1i), ...
                'magneticYY', zeros(counts - [0, 1], 'like', 1i));
            order = randperm(prod(counts));

            for chunk = {order(1:7), order(8:end)}
                parts = rooftopMoM.assembler.greensMatrix.utils.accumulateMoments(parts, moments(:, chunk{1}), ...
                    offsets(:, chunk{1}), [2, 3]);
            end

            moment = @(row) reshape(moments(row, :), counts);
            testCase.verifyEqual(parts.electric, 2 * moment(1), 'RelTol', 1e-15);
            testCase.verifyEqual(parts.magneticXX, 3 * combine_(moment(1), moment(2), moment(3), moment(4)), ...
                'RelTol', 1e-14);
            testCase.verifyEqual(parts.magneticYY, 3 * combine_(moment(1).', moment(5).', moment(6).', ...
                moment(7).').', 'RelTol', 1e-14);
        end

        function testSkipsAMissingDirection(testCase)
            counts = [6, 1];
            parts = struct('electric', zeros(counts, 'like', 1i), 'magneticXX', zeros(5, 1, 'like', 1i), ...
                'magneticYY', []);
            [moments, offsets] = deal(ones(7, 6), [0:5; zeros(1, 6)]);

            parts = rooftopMoM.assembler.greensMatrix.utils.accumulateMoments(parts, moments, offsets, [1, 1]);

            testCase.verifyEmpty(parts.magneticYY);
            testCase.verifyEqual(parts.electric, ones(6, 1));
        end
    end
end

function generator = combine_(M00, M10, M01, M11)
% combine_: The generator of the x-directed rooftops from the moments of all offsets (rows along x): at the offset d,
%   M11(d) + (M00 - M10 - M01 + M11)(d) + (M10 - M11)(d + 1) + (M01 - M11)(d - 1), with
%   (M01 - M11)(-1) = (M10 - M11)(1).

    rows = 1:(size(M00, 1) - 1);
    falling = [M10(2, :) - M11(2, :); M01(rows(1:(end - 1)), :) - M11(rows(1:(end - 1)), :)];
    generator = M11(rows, :) + (M00(rows, :) - M10(rows, :) - M01(rows, :) + M11(rows, :)) ...
        + (M10(rows + 1, :) - M11(rows + 1, :)) + falling;
end
