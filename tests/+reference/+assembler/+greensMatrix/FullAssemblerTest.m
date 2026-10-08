classdef (TestTags = {'Reference'}) FullAssemblerTest < reference.Base
% FullAssemblerTest: The generators against the reference data.

    methods (Test)
        function testGeneratorsMatchTheReference(testCase)
            for data = testCase.reference
                generators = rooftopMoM.assembler.greensMatrix.FullAssembler().generate( ...
                    testCase.createReferenceContext(data));

                % The reference is accurate to about 4e-11 on the self and neighbour entries of the magnetic part, and
                % to 1e-14 elsewhere.
                testCase.verifyEqual(generators.magneticXX(), data.magneticXX, 'RelTol', 1e-10);
                testCase.verifyEqual(generators.magneticYY(), data.magneticYY, 'RelTol', 1e-10);
                testCase.verifyEqual(generators.electric(), data.electric, 'RelTol', 1e-13);
            end
        end
    end
end
