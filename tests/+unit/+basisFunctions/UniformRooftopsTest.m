classdef UniformRooftopsTest < unit.Base
% UniformRooftopsTest: Tests of rooftopMoM.basisFunctions.UniformRooftops.

    methods (Test)
        function testCounts(testCase)
            rooftops = testCase.rooftops_(4, 3);

            testCase.verifyEqual(rooftops.countDir('x'), 9);
            testCase.verifyEqual(rooftops.countDir('y'), 8);
            testCase.verifyEqual(rooftops.count(), 17);
        end

        function testContent(testCase)
            testCase.verifyEqual(testCase.rooftops_(4, 3).content(), rooftopMoM.basisFunctions.RooftopContent.BOTH);
            testCase.verifyEqual(testCase.rooftops_(4, 1).content(), rooftopMoM.basisFunctions.RooftopContent.XONLY);
            testCase.verifyEqual(testCase.rooftops_(1, 3).content(), rooftopMoM.basisFunctions.RooftopContent.YONLY);
            testCase.verifyEqual(testCase.rooftops_(1, 1).content(), rooftopMoM.basisFunctions.RooftopContent.NONE);
        end

        function testFaceOffsetsRunWithXFastest(testCase)
            rooftops = testCase.rooftops_(4, 3);

            testCase.verifyEqual(rooftops.faceOffsets('x', [1, 3, 4, 9]), [0, 2, 0, 2; 0, 0, 1, 2]);
            testCase.verifyEqual(rooftops.faceOffsets('y', [1, 4, 5, 8]), [0, 3, 0, 3; 0, 0, 1, 1]);
        end

        function testReferenceRooftop(testCase)
            rooftops = testCase.rooftops_(4, 3, 0.1, 0.2);
            weights = rooftops.quadrature().weights2D();

            rooftopX = rooftops.referenceRooftop('x');
            testCase.verifyEqual(rooftopX.secondFaceOffset, [1; 0]);
            testCase.verifyEqual([rooftopX.length, rooftopX.height], [0.1, 0.2]);

            % A rooftop of unit height covers half of each of its two faces.
            testCase.verifyEqual(weights.' * rooftopX.amplitudes, [0.5, 0.5], 'RelTol', 1e-14);

            rooftopY = rooftops.referenceRooftop('y');
            testCase.verifyEqual(rooftopY.secondFaceOffset, [0; 1]);
            testCase.verifyEqual([rooftopY.length, rooftopY.height], [0.2, 0.1]);
        end

        function testAmplitudesRiseAndFallAlongTheDirection(testCase)
            rooftops = testCase.rooftops_(4, 3);
            [pointsU, pointsV] = rooftops.quadrature().points();

            testCase.verifyEqual(rooftops.referenceRooftop('x').amplitudes, [pointsU.', 1 - pointsU.']);
            testCase.verifyEqual(rooftops.referenceRooftop('y').amplitudes, [pointsV.', 1 - pointsV.']);
        end

        function testEdgeLengths(testCase)
            rooftops = testCase.rooftops_(4, 3, 0.1, 0.2);

            % Rooftops 1 to 9 are x-directed, 10 to 17 y-directed.
            testCase.verifyEqual(rooftops.edgeLengths([1, 9, 10, 17]), [0.2, 0.2, 0.1, 0.1]);
        end
    end

    methods (Access = private)
        function rooftops = rooftops_(testCase, countX, countY, sizeX, sizeY)
            if nargin < 4
                sizeX = 0.1;
                sizeY = 0.1;
            end

            rooftops = testCase.createRooftops(countX, countY, sizeX, sizeY, 3);
        end
    end
end
