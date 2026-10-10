classdef MeshTest < unit.Base
% MeshTest: Tests of rooftopMoM.Mesh.

    methods (Test)
        function testGetters(testCase)
            mesh = rooftopMoM.Mesh(4, 3, 0.1, 0.2);
            [countX, countY] = mesh.count();
            [sizeX, sizeY] = mesh.faceSize();

            testCase.verifyEqual([countX, countY, sizeX, sizeY], [4, 3, 0.1, 0.2]);
        end

        function testRadiusIsHalfTheDiagonal(testCase)
            % The plate is 0.4 x 0.6 m.
            testCase.verifyEqual(rooftopMoM.Mesh(4, 3, 0.1, 0.2).radius(), sqrt(0.4^2 + 0.6^2) / 2, 'RelTol', 1e-15);
        end

        function testInvalidInputsRaiseErrors(testCase)
            testCase.verifyError(@() rooftopMoM.Mesh(0, 3, 0.1, 0.1), 'MATLAB:validators:mustBePositive');
            testCase.verifyError(@() rooftopMoM.Mesh(2.5, 3, 0.1, 0.1), 'MATLAB:validators:mustBeInteger');
            testCase.verifyError(@() rooftopMoM.Mesh(4, 3, -0.1, 0.1), 'MATLAB:validators:mustBePositive');
            testCase.verifyError(@() rooftopMoM.Mesh(4, 3, 0.1, Inf), 'MATLAB:validators:mustBeFinite');
        end

        function testSingleFaceHasNoRooftop(testCase)
            testCase.verifyError(@() rooftopMoM.Mesh(1, 1, 0.1, 0.1), 'rooftopMoM:mesh:noRooftops');
        end
    end
end
