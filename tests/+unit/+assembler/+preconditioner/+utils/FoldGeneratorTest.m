classdef FoldGeneratorTest < unit.Base
% FoldGeneratorTest: Tests of rooftopMoM.assembler.preconditioner.utils.foldGenerator.

    methods (Test)
        function testFoldOfAFourByFourToeplitz(testCase)
            % The example of the header: t(-3), ..., t(3) folded with the window of T. Chan.
            generator = [-3; -2; -1; 0; 1; 2; 3] + 10;
            circulant = rooftopMoM.assembler.preconditioner.utils.foldGenerator(generator, ...
                @(offsets, count, ~) (count - abs(offsets)) / count);

            expected = [10; (3 * 11 + 7) / 4; (2 * 12 + 2 * 8) / 4; (13 + 3 * 9) / 4];
            testCase.verifyEqual(circulant, expected, 'AbsTol', 1e-14);
        end

        function testFoldsAlongBothDirections(testCase)
            % With unit weights the offset (p, q) lands at (p mod 2, q mod 3).
            generator = reshape(1:15, 3, 5);
            circulant = rooftopMoM.assembler.preconditioner.utils.foldGenerator(generator, @(offsets, ~, ~) ...
                ones(size(offsets)));

            expected = zeros(2, 3);

            for m = -1:1
                for n = -2:2
                    expected(mod(m, 2) + 1, mod(n, 3) + 1) = expected(mod(m, 2) + 1, mod(n, 3) + 1) ...
                        + generator(m + 2, n + 3);
                end
            end

            testCase.verifyEqual(circulant, expected);
        end
    end
end
