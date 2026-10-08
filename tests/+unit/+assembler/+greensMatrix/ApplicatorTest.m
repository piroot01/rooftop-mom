classdef ApplicatorTest < unit.Base
% ApplicatorTest: Tests of the three products with the Green's matrix against the dense matrix of the generators.

    properties (TestParameter)
        applicator = {'FFTApplicator', 'SplitFFTApplicator', 'LeanSplitFFTApplicator'};

        % countX, countY, sizeX, sizeY
        grid = struct( ...
            'plate', [5, 4, 0.1, 0.15], ...
            'oddPlate', [11, 7, 0.1, 0.1], ...
            'twoFacesWide', [2, 5, 0.1, 0.1], ...
            'stripX', [9, 1, 0.1, 0.02], ...
            'stripY', [1, 6, 0.03, 0.1] ...
        );
    end

    methods (Test)
        function testMatchesTheDenseMatrix(testCase, applicator, grid)
            context = testCase.createContext(grid(1), grid(2), grid(3), grid(4), 3, 2 * pi);
            rooftops = context.rooftops();
            generators = rooftopMoM.assembler.greensMatrix.FullAssembler().generate(context);
            product = feval(['rooftopMoM.assembler.greensMatrix.', applicator], rooftops, generators);

            rng(1);
            input = complex(randn(rooftops.count(), 1), randn(rooftops.count(), 1));
            expected = testCase.generateDenseMatrix(generators, rooftops) * input;

            testCase.verifyLessThan(norm(product.apply(input) - expected) / norm(expected), 1e-12);
        end

        function testMemoryIsReported(testCase, applicator)
            context = testCase.createContext(5, 4, 0.1, 0.15, 3, 2 * pi);
            generators = rooftopMoM.assembler.greensMatrix.FullAssembler().generate(context);
            product = feval(['rooftopMoM.assembler.greensMatrix.', applicator], context.rooftops(), generators);

            testCase.verifyGreaterThan(product.memory(), 0);
            testCase.verifyGreaterThan(product.workingMemory(), 16 * context.rooftops().count());
        end
    end
end
