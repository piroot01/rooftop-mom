classdef ConstantsTest < unit.Base
% ConstantsTest: Tests of rooftopMoM.assembler.Constants.

    methods (Test)
        function testConstantsAreConsistent(testCase)
            % B13: v2 mixed the pre-2019 exact mu0 with the CODATA 2018 epsilon0 and Z0.
            c = rooftopMoM.assembler.Constants.c;
            Z0 = rooftopMoM.assembler.Constants.Z0;
            mu0 = rooftopMoM.assembler.Constants.mu0;
            epsilon0 = rooftopMoM.assembler.Constants.epsilon0;

            testCase.verifyEqual(c^2 * mu0 * epsilon0, 1, 'RelTol', 4 * eps);
            testCase.verifyEqual(sqrt(mu0 / epsilon0), Z0, 'RelTol', 4 * eps);
        end

        function testMatchesCodata2018(testCase)
            % The relative standard uncertainty of mu0 and epsilon0 in CODATA 2018 is 1.5e-10.
            testCase.verifyEqual(rooftopMoM.assembler.Constants.mu0, 1.25663706212e-6, 'RelTol', 1.5e-10);
            testCase.verifyEqual(rooftopMoM.assembler.Constants.epsilon0, 8.8541878128e-12, 'RelTol', 1.5e-10);
        end
    end
end
