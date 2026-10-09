classdef GenerateOptimalOfProductTest < unit.assembler.preconditioner.Base
% GenerateOptimalOfProductTest: Tests of rooftopMoM.assembler.preconditioner.utils.generateOptimalOfProduct.

    properties (TestParameter)
        grid = {[4, 3], [5, 4], [6, 1], [1, 5]};
    end

    methods (Test)
        function testMatchesTheDenseProduct(testCase, grid)
            rng(1);
            generatorA = complex(randn(2 * grid - 1), randn(2 * grid - 1));
            generatorB = complex(randn(2 * grid - 1), randn(2 * grid - 1));

            eigenvalues = rooftopMoM.assembler.preconditioner.utils.generateOptimalOfProduct( ...
                rooftopMoM.assembler.preconditioner.utils.generatePartSpectra(generatorA), ...
                rooftopMoM.assembler.preconditioner.utils.generatePartSpectra(generatorB));

            % The optimal circulant keeps the diagonal of the matrix in the Fourier basis.
            fourier = testCase.generateFourierMatrix(grid);
            product = testCase.generateToeplitzMatrix(generatorA) * testCase.generateToeplitzMatrix(generatorB)';
            expected = reshape(diag(fourier * product * fourier'), grid);

            testCase.verifyEqual(eigenvalues, expected, 'AbsTol', 1e-12 * max(abs(expected(:))));
        end
    end
end
