classdef GeneratePartSpectraTest < unit.assembler.preconditioner.Base
% GeneratePartSpectraTest: Tests of rooftopMoM.assembler.preconditioner.utils.generatePartSpectra.

    properties (TestParameter)
        grid = {[4, 3], [6, 1], [1, 5]};
    end

    methods (Test)
        function testPartsSumToTheToeplitzMatrix(testCase, grid)
            rng(1);
            generator = complex(randn(2 * grid - 1), randn(2 * grid - 1));
            spectra = rooftopMoM.assembler.preconditioner.utils.generatePartSpectra(generator);

            % A part is F' * diag(eigenvalues) * F in the basis modulated along its skew-circulant directions.
            fourier = testCase.generateFourierMatrix(grid);
            modulationsX = {ones(grid(1), 1), exp(1i * pi * (0:(grid(1) - 1)).' / grid(1))};
            modulationsY = {ones(grid(2), 1), exp(1i * pi * (0:(grid(2) - 1)).' / grid(2))};
            matrix = zeros(prod(grid));

            for a = 1:2
                for b = 1:2
                    modulation = diag(kron(modulationsY{b}, modulationsX{a}));
                    matrix = matrix + modulation' * fourier' * diag(spectra{a, b}(:)) * fourier * modulation;
                end
            end

            testCase.verifyEqual(matrix, testCase.generateToeplitzMatrix(generator), 'AbsTol', 1e-13);
        end
    end
end
