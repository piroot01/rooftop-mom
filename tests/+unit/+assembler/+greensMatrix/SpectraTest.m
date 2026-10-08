classdef SpectraTest < unit.Base
% SpectraTest: Tests of rooftopMoM.assembler.greensMatrix.Spectra.

    methods (Test)
        function testInverseTimesSpectraIsTheIdentity(testCase)
            rng(1);
            block = @() complex(randn(4, 3), randn(4, 3));
            spectra = rooftopMoM.assembler.greensMatrix.Spectra(block(), block(), block(), block());
            inverse = spectra.invert(rooftopMoM.basisFunctions.RooftopContent.BOTH);

            % [XX, XY; YX, YY] * inverse at every frequency.
            testCase.verifyEqual(spectra.blockXX() .* inverse.blockXX() + spectra.blockXY() .* inverse.blockYX(), ...
                ones(4, 3), 'AbsTol', 1e-12);
            testCase.verifyEqual(spectra.blockXX() .* inverse.blockXY() + spectra.blockXY() .* inverse.blockYY(), ...
                zeros(4, 3), 'AbsTol', 1e-12);
            testCase.verifyEqual(spectra.blockYX() .* inverse.blockXX() + spectra.blockYY() .* inverse.blockYX(), ...
                zeros(4, 3), 'AbsTol', 1e-12);
            testCase.verifyEqual(spectra.blockYX() .* inverse.blockXY() + spectra.blockYY() .* inverse.blockYY(), ...
                ones(4, 3), 'AbsTol', 1e-12);
        end

        function testInverseOfAStripIsTheReciprocal(testCase)
            spectra = rooftopMoM.assembler.greensMatrix.Spectra([2; 4], [], [], []);
            inverse = spectra.invert(rooftopMoM.basisFunctions.RooftopContent.XONLY);

            testCase.verifyEqual(inverse.blockXX(), [0.5; 0.25]);
            testCase.verifyEmpty(inverse.blockYY());
        end

        function testGridSizeAndMemory(testCase)
            spectra = rooftopMoM.assembler.greensMatrix.Spectra([], [], [], complex(zeros(3, 5)));
            [sizeX, sizeY] = spectra.gridSize();

            testCase.verifyEqual([sizeX, sizeY], [3, 5]);
            testCase.verifyEqual(spectra.memory(), 16 * 15);
        end
    end
end
