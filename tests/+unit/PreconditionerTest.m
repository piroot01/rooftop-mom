classdef PreconditionerTest < unit.Base
% PreconditionerTest: Tests of rooftopMoM.Preconditioner.

    properties (TestParameter)
        kernel = {'strang', 'tChan', 'rChan', 'huckle'};
        builder = {'partwise', 'blockwise', 'hybrid'};
    end

    methods (Test)
        function testMatchesTheCirculantOfTheCore(testCase, kernel, builder)
            problem = testCase.problem_();
            preconditioner = rooftopMoM.Preconditioner(problem, 'kernel', kernel, 'builder', builder);

            % The default bandwidth of the huckle kernel is the whole grid of 5 x 4 faces.
            kernels = struct( ...
                'strang', @rooftopMoM.assembler.preconditioner.kernels.generateStrang, ...
                'tChan', @rooftopMoM.assembler.preconditioner.kernels.generateTChan, ...
                'rChan', @rooftopMoM.assembler.preconditioner.kernels.generateRChan, ...
                'huckle', @(generator) rooftopMoM.assembler.preconditioner.kernels.generateHuckle(generator, 5, 4));
            builders = struct( ...
                'partwise', @rooftopMoM.assembler.preconditioner.generatePartwiseSpectra, ...
                'blockwise', @rooftopMoM.assembler.preconditioner.generateBlockwiseSpectra, ...
                'hybrid', @rooftopMoM.assembler.preconditioner.generateHybridSpectra);
            expected = rooftopMoM.assembler.preconditioner.generateCirculant(problem.rooftops(), ...
                problem.generators(), kernels.(kernel), builders.(builder));

            rng(1);
            input = complex(randn(31, 1), randn(31, 1));
            testCase.verifyEqual(preconditioner.apply(input), expected.apply(input));
        end

        function testSingularApproximationRaisesError(testCase)
            testCase.verifyError(@() rooftopMoM.Preconditioner(testCase.problem_(), ...
                'modifier', @(generator) zeros(size(generator))), 'rooftopMoM:preconditioner:singular');
        end

        function testInvalidInputsRaiseErrors(testCase)
            problem = testCase.problem_();

            testCase.verifyError(@() rooftopMoM.Preconditioner(problem, 'kernel', 'jacobi'), ...
                'MATLAB:validators:mustBeMember');
            testCase.verifyError(@() rooftopMoM.Preconditioner(problem, 'kernel', 'huckle', 'bandwidth', [1, 2, 3]), ...
                'rooftopMoM:preconditioner:invalidBandwidth');
            testCase.verifyError(@() rooftopMoM.Preconditioner(problem, 'kernel', 'huckle', 'bandwidth', [0, 2]), ...
                'MATLAB:validators:mustBePositive');
        end
    end

    methods (Access = private)
        function problem = problem_(~)
            problem = rooftopMoM.Problem(rooftopMoM.Mesh(5, 4, 0.1, 0.15), 2 * pi, 'quadratureOrder', 3);
        end
    end
end
