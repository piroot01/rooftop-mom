classdef BasisFunctionsInterface < handle
    methods (Abstract)
        generate(obj);

        count = count(obj);

        indices = indices(obj);

        mesh = mesh(obj);

        memory = memory(obj);
    end
end
