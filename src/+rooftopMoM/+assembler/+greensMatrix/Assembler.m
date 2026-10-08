classdef (Abstract) Assembler
% Assembler: Interface of the assemblers of the Green's matrix generators.

    methods (Abstract)
        % generate: Returns the Generators of an AssemblerContext.
        generators = generate(obj, context);
    end
end
