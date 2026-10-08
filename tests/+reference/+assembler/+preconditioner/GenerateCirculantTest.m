classdef (TestTags = {'Reference'}) GenerateCirculantTest < reference.Base
% GenerateCirculantTest: The circulant preconditioners of v3 against v2 for every kernel and builder.

    methods (Test)
        function testPreconditionersMatchV2(testCase)
            % Rows and columns as in tests/data/generateReferenceData.m.
            kernels = {
                @rooftopMoM.assembler.preconditioner.kernels.generateStrang
                @rooftopMoM.assembler.preconditioner.kernels.generateTChan
                @rooftopMoM.assembler.preconditioner.kernels.generateRChan
            };
            builders = {
                @rooftopMoM.assembler.preconditioner.generatePartwiseSpectra
                @rooftopMoM.assembler.preconditioner.generateBlockwiseSpectra
                @rooftopMoM.assembler.preconditioner.generateHybridSpectra
            };

            for data = testCase.reference
                rooftops = testCase.createReferenceContext(data).rooftops();
                generators = rooftopMoM.assembler.greensMatrix.Generators(data.magneticXX, data.magneticYY, ...
                    data.electric, rooftops.content());

                for k = 1:numel(kernels)
                    for b = 1:numel(builders)
                        preconditioner = rooftopMoM.assembler.preconditioner.generateCirculant(rooftops, generators, ...
                            kernels{k}, builders{b});
                        expected = data.preconditioned{k, b};

                        testCase.verifyLessThan(norm(preconditioner.apply(data.vector) - expected) / norm(expected), ...
                            1e-12);
                    end
                end
            end
        end
    end
end
