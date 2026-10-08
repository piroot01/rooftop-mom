classdef GenerateMiddlePortTest < unit.Base
% GenerateMiddlePortTest: Tests of rooftopMoM.operator.utils.generateMiddlePort.

    methods (Test)
        function testStrip(testCase)
            % Six faces have the edge between the third and the fourth face in the middle; five faces take the edge
            % just before the middle, between the second and the third face.
            testCase.verifyEqual(testCase.generate_(6, 1, 'x'), 3);
            testCase.verifyEqual(testCase.generate_(5, 1, 'x'), 2);
            testCase.verifyEqual(testCase.generate_(1, 6, 'y'), 3);
        end

        function testPlate(testCase)
            % The 3x3 x-directed rooftops come first: the port is their second row. The y-directed ones form a 4x2
            % grid: the port is its first column.
            testCase.verifyEqual(testCase.generate_(4, 3, 'x'), [2, 5, 8]);
            testCase.verifyEqual(testCase.generate_(4, 3, 'y'), 9 + (1:4));
        end
    end

    methods (Access = private)
        function port = generate_(testCase, countX, countY, direction)
            rooftops = testCase.createRooftops(countX, countY, 0.1, 0.1, 3);
            port = rooftopMoM.operator.utils.generateMiddlePort(rooftops, direction);
        end
    end
end
