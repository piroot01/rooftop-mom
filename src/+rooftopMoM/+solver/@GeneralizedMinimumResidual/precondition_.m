function output = precondition_(obj, input)
% precondition_: Returns the preconditioner applied to the input, the input itself without one.

    output = input;

    if ~isempty(obj.preconditioner_)
        output = obj.preconditioner_(input);
    end
end
