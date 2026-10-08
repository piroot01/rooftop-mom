classdef Solver
% Solver: Iterative solver of the system of a Problem, or of any other linear system.
%       solver = rooftopMoM.Solver('method', 'cg', 'tolerance', 1e-6);
%       result = solver.solve(operator, excitation, 'preconditioner', preconditioner);
%   The operator and the preconditioner may also be a matrix, a function handle or any object with the method apply,
%   and the right-hand side a plain vector. Outside of this package the solver needs only rooftopMoM.Solver,
%   rooftopMoM.Result and the packages rooftopMoM.solver and rooftopMoM.logging.

    properties (Access = private)
        method_ char;
        tolerance_ double;
        maxIterationCount_ double;
        restart_ double;
    end

    methods
        function obj = Solver(options)
        % Solver: Creates the solver.
        %
        %   Optional name-value inputs:
        %       'method': 'cg', the conjugate gradient method for complex symmetric matrices like the Green's matrix
        %           (COCG), or 'gmres', the restarted GMRES of MATLAB for any matrix (default 'cg').
        %       'tolerance': Relative residual norm to stop at, norm(b - A * x) / norm(b) (default 1e-6).
        %       'maxIterationCount': Maximum number of iterations (default 1000).
        %       'restart': Iterations between the restarts of 'gmres', Inf for none (default 50).

            arguments
                options.method (1, :) char {mustBeMember(options.method, {'cg', 'gmres'})} = 'cg'
                options.tolerance (1, 1) {mustBeReal, mustBePositive, mustBeLessThan(options.tolerance, 1)} = 1e-6
                options.maxIterationCount (1, 1) {mustBeInteger, mustBePositive} = 1000
                options.restart (1, 1) {mustBePositive, mustBeIntegerOrInf_} = 50
            end

            obj.method_ = options.method;
            obj.tolerance_ = double(options.tolerance);
            obj.maxIterationCount_ = double(options.maxIterationCount);
            obj.restart_ = double(options.restart);
        end
    end

    methods
        function result = solve(obj, operator, rightHandSide, options)
        % solve: Solves operator * x = rightHandSide.
        %
        %   Inputs:
        %       operator: rooftopMoM.Operator, matrix, function handle or object with the method apply.
        %       rightHandSide: rooftopMoM.Excitation or numeric vector.
        %
        %   Optional name-value inputs:
        %       'preconditioner': rooftopMoM.Preconditioner, matrix, function handle or object with the method apply,
        %           an approximate inverse of the operator (default none).
        %       'initialGuess': Start vector (default zeros).
        %
        %   Output:
        %       result: rooftopMoM.Result.

            arguments
                obj
                operator {mustBeOperator_}
                rightHandSide {mustBeRightHandSide_}
                options.preconditioner {mustBeOperatorOrEmpty_} = []
                options.initialGuess {mustBeNumeric, mustBeFinite} = []
            end

            if isa(rightHandSide, 'rooftopMoM.Excitation')
                excitation = rightHandSide;
                vector = rightHandSide.vector();
            else
                excitation = [];
                vector = double(rightHandSide(:));
            end

            count = numel(vector);
            checkSize_(operator, count, 'operator');
            checkSize_(options.preconditioner, count, 'preconditioner');

            if ~isempty(options.initialGuess) && numel(options.initialGuess) ~= count
                error('rooftopMoM:solver:sizeMismatch', 'The initial guess must have %d elements.', count);
            end

            method = obj.createMethod_(unwrap_(options.preconditioner));
            output = method.solve(unwrap_(operator), vector, double(options.initialGuess(:)));
            result = rooftopMoM.Result(output, excitation);
        end
    end

    methods (Access = private)
        function method = createMethod_(obj, preconditioner)
        % createMethod_: Returns the rooftopMoM.solver.IterativeMethod of the chosen method.

            switch obj.method_
                case 'cg'
                    if isempty(preconditioner)
                        method = rooftopMoM.solver.ConjugateGradient(obj.tolerance_, obj.maxIterationCount_);
                    else
                        method = rooftopMoM.solver.PreconditionedConjugateGradient(preconditioner, obj.tolerance_, ...
                            obj.maxIterationCount_);
                    end
                case 'gmres'
                    method = rooftopMoM.solver.GeneralizedMinimumResidual(obj.tolerance_, obj.maxIterationCount_, ...
                        obj.restart_, preconditioner);
            end
        end
    end
end

function operator = unwrap_(operator)
    % The iterations use the core product, which skips the checks of the API.
    if isa(operator, 'rooftopMoM.Product')
        operator = operator.applicator();
    end
end

function checkSize_(operator, count, name)
    if isnumeric(operator) && ~isempty(operator)
        matches = isequal(size(operator), [count, count]);
    elseif isa(operator, 'rooftopMoM.Product')
        matches = operator.problem().rooftopCount() == count;
    else
        % A function handle or another object reveals its size only when applied.
        matches = true;
    end

    if ~matches
        error('rooftopMoM:solver:sizeMismatch', 'The %s does not match the %d unknowns of the right-hand side.', ...
            name, count);
    end
end

function mustBeOperator_(operator)
    isOperator = isnumeric(operator) || isa(operator, 'function_handle') || ...
        (isstruct(operator) && isfield(operator, 'apply')) || (isobject(operator) && ismethod(operator, 'apply'));

    if ~isOperator || isempty(operator)
        error('rooftopMoM:solver:invalidOperator', ['An operator is a matrix, a function handle or an object with ' ...
            'the method apply.']);
    end
end

function mustBeOperatorOrEmpty_(operator)
    if ~isempty(operator)
        mustBeOperator_(operator);
    end
end

function mustBeRightHandSide_(rightHandSide)
    isVector = isnumeric(rightHandSide) && isvector(rightHandSide) && allfinite(rightHandSide);

    if ~isVector && ~isa(rightHandSide, 'rooftopMoM.Excitation')
        error('rooftopMoM:solver:invalidRightHandSide', ['The right-hand side is a rooftopMoM.Excitation or a ' ...
            'numeric vector of finite values.']);
    end
end

function mustBeIntegerOrInf_(value)
    if ~isinf(value) && value ~= round(value)
        error('rooftopMoM:solver:invalidRestart', 'The restart is a positive integer or Inf.');
    end
end
