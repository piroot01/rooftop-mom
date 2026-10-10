classdef OperatorTest < unit.Base
% OperatorTest: Tests of rooftopMoM.Operator and of the checks of rooftopMoM.Product.

    properties (TestParameter)
        product = {'fft', 'splitFFT', 'leanSplitFFT'};
    end

    methods (Test)
        function testMatchesTheDenseMatrix(testCase, product)
            problem = testCase.problem_();
            operator = rooftopMoM.Operator(problem, 'product', product);
            matrix = testCase.generateDenseMatrix(problem.generators(), problem.rooftops());
            input = testCase.input_(problem);

            testCase.verifyLessThan(norm(operator.apply(input) - matrix * input) / norm(matrix * input), 1e-12);
        end

        function testModifierSelectsAPartOfTheMatrix(testCase)
            problem = testCase.problem_();
            operator = rooftopMoM.Operator(problem, 'modifier', @imag);
            matrix = imag(testCase.generateDenseMatrix(problem.generators(), problem.rooftops()));
            input = testCase.input_(problem);

            testCase.verifyLessThan(norm(operator.apply(input) - matrix * input) / norm(matrix * input), 1e-12);
        end

        function testApplyChecksTheInput(testCase)
            operator = rooftopMoM.Operator(testCase.problem_());

            testCase.verifyError(@() operator.apply(ones(3, 1)), 'rooftopMoM:product:sizeMismatch');
            testCase.verifyError(@() operator.apply(ones(1, 31)), 'rooftopMoM:product:sizeMismatch');
            testCase.verifyError(@() operator.apply(repmat('a', 31, 1)), 'rooftopMoM:product:sizeMismatch');
        end

        function testUnknownProductRaisesError(testCase)
            testCase.verifyError(@() rooftopMoM.Operator(testCase.problem_(), 'product', 'dense'), ...
                'MATLAB:validators:mustBeMember');
        end
    end

    methods (Access = private)
        function problem = problem_(~)
            % At the wavenumber 4 the faces of 0.15 m stay below a tenth of the wavelength.
            problem = rooftopMoM.Problem(rooftopMoM.Mesh(5, 4, 0.1, 0.15), 4, 'quadratureOrder', 3);
        end

        function input = input_(~, problem)
            rng(1);
            input = complex(randn(problem.rooftopCount(), 1), randn(problem.rooftopCount(), 1));
        end
    end
end
