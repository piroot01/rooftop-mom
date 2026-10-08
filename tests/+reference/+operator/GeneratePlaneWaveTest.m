classdef (TestTags = {'Reference'}) GeneratePlaneWaveTest < reference.Base
% GeneratePlaneWaveTest: The plane wave excitation against the reference data.

    methods (Test)
        function testExcitationMatchesTheReference(testCase)
            for data = testCase.reference
                rooftops = testCase.createReferenceContext(data).rooftops();
                excitation = rooftopMoM.operator.generatePlaneWave(rooftops, data.wavenumber, data.direction, ...
                    data.polarization, data.intensity);

                testCase.verifyLessThan(norm(excitation - data.planeWave) / norm(data.planeWave), 1e-13);
            end
        end
    end
end
