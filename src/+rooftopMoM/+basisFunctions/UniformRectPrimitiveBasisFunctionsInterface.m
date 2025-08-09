classdef UniformRectPrimitiveBasisFunctionsInterface < handle
    methods (Abstract)
        generate(obj);

        count = count(obj, dir);

        indices = indices(obj);

        primitive = primitive(obj);

        memory = memory(obj);
    end
end
