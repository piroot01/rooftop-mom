classdef BasisFunctionsInterface < handle
    methods (Abstract)
        generate(obj, genome);

        count = count(obj);

        indices = indices(obj);

        mesh = mesh(obj);
    end
end
