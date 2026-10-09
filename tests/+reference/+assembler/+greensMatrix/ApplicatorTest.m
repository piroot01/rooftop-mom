classdef (TestTags = {'Reference'}) ApplicatorTest < reference.Base
% ApplicatorTest: The products against the reference data.

    methods (Test)
        function testProductsMatchTheReference(testCase)
            names = {'FFTApplicator', 'SplitFFTApplicator', 'LeanSplitFFTApplicator'};
            fields = {'productFFT', 'productSplit', 'productLean'};

            for data = testCase.reference
                % The reference generators as input, so only the products are compared.
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

        function testProductMatchesTheReferenceDenseMatrix(testCase)
            % The reference dense matrix is assembled without generators, with the quadrature of the reference
            % generators, and is accurate to about 4e-11.
            data = testCase.reference(1);
            rooftops = testCase.createReferenceContext(data).rooftops();
            generators = rooftopMoM.assembler.greensMatrix.Generators(data.magneticXX, data.magneticYY, ...
                data.electric, rooftops.content());
            product = testCase.createProduct('FFTApplicator', rooftops, generators);
            expected = data.impedanceMatrix * data.vector;

            testCase.verifyLessThan(norm(product.apply(data.vector) - expected) / norm(expected), 1e-10);
        end
    end
end
