classdef FullAssemblerTest < unit.Base
% FullAssemblerTest: Tests of rooftopMoM.assembler.greensMatrix.FullAssembler.

    methods (Test)
        function testGeneratorSizes(testCase)
            generators = testCase.generate_(5, 4);

            testCase.verifySize(generators.magneticXX(), [4, 4]);
            testCase.verifySize(generators.magneticYY(), [5, 3]);
            testCase.verifySize(generators.electric(), [5, 4]);
            testCase.verifyEqual(generators.content(), rooftopMoM.basisFunctions.RooftopContent.BOTH);
        end

        function testStripHasOneDirection(testCase)
            generators = testCase.generate_(6, 1);

            testCase.verifySize(generators.magneticXX(), [5, 1]);
            testCase.verifyEmpty(generators.magneticYY());
            testCase.verifyEqual(generators.content(), rooftopMoM.basisFunctions.RooftopContent.XONLY);
        end

        function testEntriesDoNotDependOnTheGridSize(testCase)
            % B16: on grids two faces wide, v2 added the self-face integral to the rooftop one row up.
            narrow = testCase.generate_(2, 3);
            wide = testCase.generate_(5, 3);
            wideXX = wide.magneticXX();
            wideYY = wide.magneticYY();
            wideElectric = wide.electric();

            testCase.verifyEqual(narrow.magneticXX(), wideXX(1, :), 'RelTol', 1e-13);
            testCase.verifyEqual(narrow.magneticYY(), wideYY(1:2, :), 'RelTol', 1e-13);
            testCase.verifyEqual(narrow.electric(), wideElectric(1:2, :), 'RelTol', 1e-13);
        end

        function testRealPartIsSymmetricInTheFaceSizes(testCase)
            % Swapping the face sizes maps the x-directed rooftops onto the y-directed ones.
            generators = testCase.generate_(5, 4, 0.1, 0.15);
            swapped = testCase.generate_(4, 5, 0.15, 0.1);

            testCase.verifyEqual(swapped.magneticYY(), generators.magneticXX().', 'RelTol', 1e-13);
            testCase.verifyEqual(swapped.electric(), generators.electric().', 'RelTol', 1e-13);
        end

        function testInfoMessageIsLogged(testCase)
            sink = rooftopMoM.logging.sinks.MemorySink();
            rooftopMoM.logging.Logger.instance().configure(rooftopMoM.logging.LogLevel.INFO, {sink});

            testCase.generate_(3, 3);

            lines = sink.lines();
            testCase.verifyNumElements(lines, 1);
            testCase.verifySubstring(lines{1}, '[FullAssembler.generate]: Assembled the generators of 12 rooftops');
        end
    end

    methods (Access = private)
        function generators = generate_(testCase, countX, countY, sizeX, sizeY)
            if nargin < 4
                sizeX = 0.1;
                sizeY = 0.1;
            end

            context = testCase.createContext(countX, countY, sizeX, sizeY, 3, 2 * pi);
            generators = rooftopMoM.assembler.greensMatrix.FullAssembler().generate(context);
        end
    end
end
