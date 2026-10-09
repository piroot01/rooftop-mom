classdef (TestTags = {'Reference'}) GenerateCirculantTest < reference.Base
% GenerateCirculantTest: The partwise circulant preconditioners against the reference data.

    methods (Test)
        function testPreconditionersMatchTheReference(testCase)
            % In the order of tests/data/generateReferenceData.m.
            kernels = {
                @rooftopMoM.assembler.preconditioner.kernels.generateStrang
                @rooftopMoM.assembler.preconditioner.kernels.generateTChan
            };

            for data = testCase.reference
                rooftops = testCase.createReferenceContext(data).rooftops();
                generators = rooftopMoM.assembler.greensMatrix.Generators(data.magneticXX, data.magneticYY, ...
                    data.electric, rooftops.content());

                for k = 1:numel(kernels)
                    preconditioner = rooftopMoM.assembler.preconditioner.generateCirculant(rooftops, generators, ...
                        kernels{k}, @rooftopMoM.assembler.preconditioner.generatePartwiseSpectra);
                    expected = data.preconditioned{k};

                    testCase.verifyLessThan(norm(preconditioner.apply(data.vector) - expected) / norm(expected), 1e-12);
                end
            end
        end
    end
end
