classdef (TestTags = {'Reference'}) FullAssemblerTest < reference.Base
% FullAssemblerTest: The generators of v3 against v2 on the reference grids.

    methods (Test)
        function testGeneratorsMatchV2(testCase)
            for data = testCase.reference
                generators = rooftopMoM.assembler.greensMatrix.FullAssembler().generate( ...
                    testCase.createReferenceContext(data));

                % B10: v2 computed the distances as sqrt(|a|^2 + |b|^2 - 2 a.b). Coincident points of the shifted
                % rooftop faces got R near 1e-9 instead of 0, which cost up to 4e-11 on the self and neighbour
                % entries of the magnetic part. All other entries agree to 1e-14.
                testCase.verifyEqual(generators.magneticXX(), data.magneticXX, 'RelTol', 1e-10);
                testCase.verifyEqual(generators.magneticYY(), data.magneticYY, 'RelTol', 1e-10);
                testCase.verifyEqual(generators.electric(), data.electric, 'RelTol', 1e-13);
            end
        end
    end
end
