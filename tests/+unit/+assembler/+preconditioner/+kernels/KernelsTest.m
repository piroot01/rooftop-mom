classdef KernelsTest < unit.assembler.preconditioner.Base
% KernelsTest: Tests of the circulant approximations in rooftopMoM.assembler.preconditioner.kernels.

    properties (TestParameter)
        kernel = struct( ...
            'strang', @rooftopMoM.assembler.preconditioner.kernels.generateStrang, ...
            'tChan', @rooftopMoM.assembler.preconditioner.kernels.generateTChan, ...
            'rChan', @rooftopMoM.assembler.preconditioner.kernels.generateRChan, ...
            'huckle', @(generator) rooftopMoM.assembler.preconditioner.kernels.generateHuckle(generator, 3, 2));
    end

    methods (Test)
        function testReversedGeneratorGivesTheTransposedCirculant(testCase, kernel)
            % g(-p, -q) in place of g(p, q) transposes the matrix, so the circulant must become c(-k, -l). This keeps
            % symmetric matrices and transposed coupling blocks symmetric.
            rng(1);
            generator = randn(11, 9);
            circulant = kernel(generator);

            testCase.verifyEqual(kernel(rot90(generator, 2)), circulant([1, 6:(-1):2], [1, 5:(-1):2]), ...
                'AbsTol', 1e-14);
        end

        function testStrangAndTChanKeepACirculant(testCase)
            % The Toeplitz matrix of a circulant generator is the circulant, its own best approximation.
            rng(1);
            circulant = randn(6, 5);
            generator = circulant(mod(-5:5, 6) + 1, mod(-4:4, 5) + 1);

            testCase.verifyEqual(rooftopMoM.assembler.preconditioner.kernels.generateStrang(generator), circulant, ...
                'AbsTol', 1e-14);
            testCase.verifyEqual(rooftopMoM.assembler.preconditioner.kernels.generateTChan(generator), circulant, ...
                'AbsTol', 1e-14);
        end

        function testStrangCopiesTheCentralOffsets(testCase)
            generator = reshape(1:20, 5, 4);
            circulant = rooftopMoM.assembler.preconditioner.kernels.generateStrang( ...
                rooftopMoM.assembler.preconditioner.utils.extendEvenly(generator));

            testCase.verifyEqual(circulant(1:3, 1:3), generator(1:3, 1:3));
            testCase.verifyEqual(circulant(4:5, 1), generator([3, 2], 1));
            testCase.verifyEqual(circulant(1, 4), generator(1, 2));
        end

        function testRChanAddsTheReflectionInEveryDirection(testCase)
            generator = reshape(1:12, 4, 3);
            circulant = rooftopMoM.assembler.preconditioner.kernels.generateRChan( ...
                rooftopMoM.assembler.preconditioner.utils.extendEvenly(generator));

            % Offsets (i, j) counted from 0: (1, 0) gains (3, 0), (0, 1) gains (0, 2), (1, 1) gains all three.
            testCase.verifyEqual(circulant(1, 1), generator(1, 1));
            testCase.verifyEqual(circulant(2, 1), generator(2, 1) + generator(4, 1));
            testCase.verifyEqual(circulant(1, 2), generator(1, 2) + generator(1, 3));
            testCase.verifyEqual(circulant(2, 2), generator(2, 2) + generator(4, 2) + generator(2, 3) + generator(4, 3));
        end

        function testHuckleKeepsTheZeroOffset(testCase)
            generator = reshape(1:20, 5, 4);
            circulant = rooftopMoM.assembler.preconditioner.kernels.generateHuckle( ...
                rooftopMoM.assembler.preconditioner.utils.extendEvenly(generator), 1, 1);

            expected = zeros(5, 4);
            expected(1, 1) = generator(1, 1);
            testCase.verifyEqual(circulant, expected);
        end

        function testHuckleDampsAndWraps(testCase)
            generator = reshape(1:20, 5, 4);
            circulant = rooftopMoM.assembler.preconditioner.kernels.generateHuckle( ...
                rooftopMoM.assembler.preconditioner.utils.extendEvenly(generator), 2, 1);

            % Offset 1 in x, damped by 1 - 1 / 2, sits at offsets 1 and 4.
            testCase.verifyEqual(circulant([2, 5], 1), [0.5; 0.5] * generator(2, 1));
            testCase.verifyEqual(circulant([3, 4], 1), [0; 0]);
        end

        function testEmptyGeneratorGivesEmptyCirculant(testCase)
            testCase.verifyEmpty(rooftopMoM.assembler.preconditioner.kernels.generateStrang([]));
            testCase.verifyEmpty(rooftopMoM.assembler.preconditioner.kernels.generateTChan([]));
            testCase.verifyEmpty(rooftopMoM.assembler.preconditioner.kernels.generateRChan([]));
            testCase.verifyEmpty(rooftopMoM.assembler.preconditioner.kernels.generateHuckle([], 2, 2));
            testCase.verifyEmpty(rooftopMoM.assembler.preconditioner.kernels.generateSuperoptimal([]));
        end

        function testSuperoptimalMatchesItsDefinition(testCase)
            % C = c(T * T') / c(T') with the eigenvalues of c, the diagonal of the matrix in the Fourier basis.
            rng(1);
            generator = complex(randn(7, 5), randn(7, 5));
            matrix = testCase.generateToeplitzMatrix(generator);
            fourier = testCase.generateFourierMatrix([4, 3]);
            expected = diag(fourier * (matrix * matrix') * fourier') ./ conj(diag(fourier * matrix * fourier'));

            circulant = rooftopMoM.assembler.preconditioner.kernels.generateSuperoptimal(generator);

            testCase.verifyEqual(reshape(fft2(circulant), [], 1), expected, 'AbsTol', 1e-12 * max(abs(expected)));
        end

        function testSuperoptimalOfAnEvenGeneratorIsSymmetric(testCase)
            rng(1);
            generator = rooftopMoM.assembler.preconditioner.utils.extendEvenly(complex(randn(5, 4), randn(5, 4)));
            circulant = rooftopMoM.assembler.preconditioner.kernels.generateSuperoptimal(generator);

            testCase.verifyEqual(circulant, circulant([1, 5:(-1):2], [1, 4:(-1):2]), ...
                'AbsTol', 1e-13 * max(abs(circulant(:))));
        end
    end
end
