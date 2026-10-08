classdef (TestTags = {'Reference'}) ApplicatorTest < reference.Base
% ApplicatorTest: The products of v3 against v2 on the reference grids.

    methods (Test)
        function testProductsMatchV2(testCase)
            names = {'FFTApplicator', 'SplitFFTApplicator', 'LeanSplitFFTApplicator'};
            fields = {'productFFT', 'productSplit', 'productLean'};

            for data = testCase.reference
                % The v2 generators as input, so only the products are compared.
                rooftops = testCase.createReferenceContext(data).rooftops();
                generators = rooftopMoM.assembler.greensMatrix.Generators(data.magneticXX, data.magneticYY, ...
                    data.electric, rooftops.content());

                for n = 1:numel(names)
                    product = testCase.createProduct(names{n}, rooftops, generators);
                    expected = data.(fields{n});
                    testCase.verifyLessThan(norm(product.apply(data.vector) - expected) / norm(expected), 1e-13);
                end
            end
        end

        function testWholeChainMatchesTheDenseMatrixOfV2(testCase)
            % The dense matrix of v2 comes from its separate impedance matrix assembler.
            data = testCase.reference(1);
            context = testCase.createReferenceContext(data);
            generators = rooftopMoM.assembler.greensMatrix.FullAssembler().generate(context);
            product = testCase.createProduct('FFTApplicator', context.rooftops(), generators);
            expected = data.impedanceMatrix * data.vector;

            % B10 limits v2 to about 4e-11 (its dense and fast paths agree to 3.6e-11).
            testCase.verifyLessThan(norm(product.apply(data.vector) - expected) / norm(expected), 1e-10);
        end
    end
end
