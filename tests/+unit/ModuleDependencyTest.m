classdef ModuleDependencyTest < unit.Base
% ModuleDependencyTest: The modules of the package use each other in one direction only and need base MATLAB only, so
%   the solver or the assembler can be copied into another application with just the modules below them.

    properties (TestParameter)
        % Each module with the modules it may use.
        module = struct( ...
            'logging', struct('name', 'logging', 'allowed', {{}}), ...
            'solver', struct('name', 'solver', 'allowed', {{'logging'}}), ...
            'mesh', struct('name', 'mesh', 'allowed', {{}}), ...
            'basisFunctions', struct('name', 'basisFunctions', 'allowed', {{'mesh', 'logging'}}), ...
            'assembler', struct('name', 'assembler', 'allowed', {{'basisFunctions', 'mesh', 'logging'}}), ...
            'operator', struct('name', 'operator', 'allowed', {{'assembler', 'basisFunctions', 'mesh', 'logging'}}), ...
            'resources', struct('name', 'resources', 'allowed', {{'logging'}}));
    end

    methods (Test)
        function testUsesOnlyTheAllowedModules(testCase, module)
            root = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'src', '+rooftopMoM');
            files = dir(fullfile(root, ['+', module.name], '**', '*.m'));
            [required, products] = matlab.codetools.requiredFilesAndProducts(fullfile({files.folder}, {files.name}));

            used = unique(cellfun(@moduleOf_, extractAfter(required, [root, filesep]), 'UniformOutput', false));

            testCase.verifyEmpty(setdiff(used, [{module.name}, module.allowed]), ...
                sprintf('The module %s uses modules it may not.', module.name));
            testCase.verifyEqual({products.Name}, {'MATLAB'});
        end
    end
end

function name = moduleOf_(path)
    % The first folder below the package is the module; the files directly in the package are the API.
    parts = strsplit(path, {'/', '\'});

    if startsWith(parts{1}, '+')
        name = parts{1}(2:end);
    else
        name = 'API';
    end
end
