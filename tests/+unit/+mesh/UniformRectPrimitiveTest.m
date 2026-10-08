classdef UniformRectPrimitiveTest < unit.Base
% UniformRectPrimitiveTest: Tests of rooftopMoM.mesh.UniformRectPrimitive.

    methods (Test)
        function testGetters(testCase)
            primitive = rooftopMoM.mesh.UniformRectPrimitive(4, 3, 0.1, 0.2);

            testCase.verifyEqual([primitive.countX(), primitive.countY()], [4, 3]);
            testCase.verifyEqual([primitive.sizeX(), primitive.sizeY()], [0.1, 0.2]);
            testCase.verifyEqual(primitive.count(), 12);
        end
    end
end
