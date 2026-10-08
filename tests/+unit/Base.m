classdef (Abstract) Base < matlab.unittest.TestCase
% Base: Common fixture of the unit tests.
%   It puts src on the path, so a test class also runs with only the folder tests on the path. The shared logger is
%   reset before and after every test, so no test depends on the logging state left by another one.

    methods (TestClassSetup)
        function addSourcePath(testCase)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'src')));
        end
    end

    methods (TestMethodSetup)
        function resetLogger(testCase)
            logger = rooftopMoM.logging.Logger.instance();
            logger.reset();
            testCase.addTeardown(@() logger.reset());
        end
    end

    methods (Access = protected)
        function context = createContext(~, countX, countY, sizeX, sizeY, order, wavenumber, memoryBytes)
        % createContext: Returns the AssemblerContext of a grid, by default with a 64 MB chunk budget.

            if nargin < 8
                memoryBytes = 2^26;
            end

            primitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
            quadrature = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);
            rooftops = rooftopMoM.basisFunctions.UniformRooftops(primitive, quadrature);
            context = rooftopMoM.assembler.AssemblerContext(rooftops, wavenumber, memoryBytes);
        end
    end
end
