classdef KernelsTest < unit.Base
% KernelsTest: Tests of the circulant approximations in rooftopMoM.assembler.preconditioner.kernels.

    methods (Test)
        function testStrangAndTChanKeepACirculant(testCase)
            % A symmetric circulant generator t(k) = t(n - k) is its own best circulant approximation.
            circulant = testCase.symmetricCirculant_(6, 5);

            testCase.verifyEqual(rooftopMoM.assembler.preconditioner.kernels.generateStrang(circulant), circulant, ...
                'AbsTol', 1e-14);
            testCase.verifyEqual(rooftopMoM.assembler.preconditioner.kernels.generateTChan(circulant), circulant, ...
                'AbsTol', 1e-14);
        end

        function testStrangCopiesTheCentralOffsets(testCase)
            generator = reshape(1:20, 5, 4);
            circulant = rooftopMoM.assembler.preconditioner.kernels.generateStrang(generator);

            testCase.verifyEqual(circulant(1:3, 1:3), generator(1:3, 1:3));
            testCase.verifyEqual(circulant(4:5, 1), generator([3, 2], 1));
            testCase.verifyEqual(circulant(1, 4), generator(1, 2));
        end

        function testRChanAddsThePointReflection(testCase)
            generator = reshape(1:12, 4, 3);
            circulant = rooftopMoM.assembler.preconditioner.kernels.generateRChan(generator);

            testCase.verifyEqual(circulant(1, :), generator(1, :));
            testCase.verifyEqual(circulant(:, 1), generator(:, 1));
            testCase.verifyEqual(circulant(2, 2), generator(2, 2) + generator(4, 3));
        end

        function testHuckleKeepsTheZeroOffset(testCase)
            % B17: v2 never set the zero offset, so its circulant missed the diagonal.
            generator = reshape(1:20, 5, 4);
            circulant = rooftopMoM.assembler.preconditioner.kernels.generateHuckle(generator, 1, 1);

            expected = zeros(5, 4);
            expected(1, 1) = generator(1, 1);
            testCase.verifyEqual(circulant, expected);
        end

        function testHuckleDampsAndWraps(testCase)
            generator = reshape(1:20, 5, 4);
            circulant = rooftopMoM.assembler.preconditioner.kernels.generateHuckle(generator, 2, 1);

            % Offset 1 in x, damped by 1 - 1 / 2, sits at offsets 1 and 4.
            testCase.verifyEqual(circulant([2, 5], 1), [0.5; 0.5] * generator(2, 1));
            testCase.verifyEqual(circulant([3, 4], 1), [0; 0]);
        end

        function testEmptyGeneratorGivesEmptyCirculant(testCase)
            testCase.verifyEmpty(rooftopMoM.assembler.preconditioner.kernels.generateStrang([]));
            testCase.verifyEmpty(rooftopMoM.assembler.preconditioner.kernels.generateTChan([]));
            testCase.verifyEmpty(rooftopMoM.assembler.preconditioner.kernels.generateRChan([]));
            testCase.verifyEmpty(rooftopMoM.assembler.preconditioner.kernels.generateHuckle([], 2, 2));
        end
    end

    methods (Access = private)
        function generator = symmetricCirculant_(~, countX, countY)
            rng(1);
            values = randn(countX, countY);
            p = mod(-(0:(countX - 1)), countX) + 1;
            q = mod(-(0:(countY - 1)), countY) + 1;

            % Averaging with the wrapped offsets makes t(k) = t(n - k) in both directions.
            generator = (values + values(p, :) + values(:, q) + values(p, q)) / 4;
        end
    end
end
