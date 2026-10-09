classdef ExtendEvenlyTest < unit.Base
% ExtendEvenlyTest: Tests of rooftopMoM.assembler.preconditioner.utils.extendEvenly.

    methods (Test)
        function testMirrorsBothOffsets(testCase)
            % Rows and columns hold the offsets -1, 0, 1, and g(p, q) = g(|p|, |q|).
            twoSided = rooftopMoM.assembler.preconditioner.utils.extendEvenly([1, 2; 3, 4]);

            testCase.verifyEqual(twoSided, [4, 3, 4; 2, 1, 2; 4, 3, 4]);
        end

        function testEmptyStaysEmpty(testCase)
            testCase.verifyEmpty(rooftopMoM.assembler.preconditioner.utils.extendEvenly([]));
        end
    end
end
