classdef (Abstract) Base < unit.Base
% Base: Common fixture of the reference tests. They compare the package with the reference data in
%   tests/data/reference.mat, computed by tests/data/generateReferenceData.m.

    properties (SetAccess = private)
        % 1xC struct array, one element per reference grid.
        reference struct;
    end

    methods (TestClassSetup)
        function loadReference(testCase)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            data = load(fullfile(root, 'tests', 'data', 'reference.mat'), 'reference');
            testCase.reference = data.reference;
        end
    end

    methods (Access = protected)
        function context = createReferenceContext(testCase, data)
        % createReferenceContext: Returns the AssemblerContext of one reference grid.

            context = testCase.createContext(data.countX, data.countY, data.sizeX, data.sizeY, data.order, ...
                data.wavenumber);
        end
    end
end
