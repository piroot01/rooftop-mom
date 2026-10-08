function handle = convertToHandle(operator)
% convertToHandle: Returns a function handle that applies an operator.
%
%   Input:
%       operator: Matrix, function handle, or any object (or struct) with the method (or field) apply.
%
%   Output:
%       handle: Function handle x -> operator applied to x.

    if isa(operator, 'function_handle')
        handle = operator;
    elseif isnumeric(operator)
        handle = @(input) operator * input;
    else
        handle = @(input) operator.apply(input);
    end
end
