function plan = buildfile
% buildfile: Build plan of the rooftopMoM package, run as buildtool <task> from the repository root.
%
%   Output:
%       plan: matlab.buildtool.Plan with the tasks:
%           check: fails on code issues of the severity error in src and tests;
%           test: runs the tests without the tag Slow (default);
%           testAll: runs all tests;
%           coverage: runs the tests without the tag Slow and writes the HTML coverage report of src into
%               build/coverage.

    plan = buildplan(localfunctions);
    plan.DefaultTasks = ["check", "test"];
end

function checkTask(~)
% checkTask: Fails on code issues of the severity error in src and tests.

    root = fileparts(mfilename('fullpath'));
    issues = codeIssues(fullfile(root, ["src", "tests"])).Issues;
    errors = issues(issues.Severity == "error", :);

    if ~isempty(errors)
        disp(errors);
        error('rooftopMoM:build:codeIssues', 'There are %d code issues of the severity error.', height(errors));
    end
end

function testTask(~)
% testTask: Runs the tests without the tag Slow.

    runTests_(~matlab.unittest.selectors.HasTag('Slow'), []);
end

function testAllTask(~)
% testAllTask: Runs all tests.

    runTests_([], []);
end

function coverageTask(~)
% coverageTask: Runs the tests without the tag Slow and writes the HTML coverage report of src into build/coverage.

    root = fileparts(mfilename('fullpath'));
    report = matlab.unittest.plugins.codecoverage.CoverageReport(fullfile(root, 'build', 'coverage'));
    plugin = matlab.unittest.plugins.CodeCoveragePlugin.forFolder(fullfile(root, 'src'), ...
        'IncludingSubfolders', true, 'Producing', report);
    runTests_(~matlab.unittest.selectors.HasTag('Slow'), plugin);
end

function runTests_(selector, plugin)
% runTests_: Runs the tests of the folder tests and fails on any failure.
%
%   Inputs:
%       selector: matlab.unittest.selectors selector of the tests to run, [] for all tests.
%       plugin: matlab.unittest.plugins plugin added to the runner, [] for none.

    root = fileparts(mfilename('fullpath'));
    addpath(fullfile(root, 'src'), fullfile(root, 'tests'));

    suite = matlab.unittest.TestSuite.fromFolder(fullfile(root, 'tests'), 'IncludingSubfolders', true);

    if ~isempty(selector)
        suite = suite.selectIf(selector);
    end

    runner = matlab.unittest.TestRunner.withTextOutput();

    if ~isempty(plugin)
        runner.addPlugin(plugin);
    end

    results = runner.run(suite);
    disp(results);
    assertSuccess(results);
end
