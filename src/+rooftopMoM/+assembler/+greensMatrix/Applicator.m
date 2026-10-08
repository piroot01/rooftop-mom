classdef (Abstract) Applicator
% Applicator: Base of the matrix-free products with the Green's matrix of a uniform grid. The rooftop coefficients are
%   ordered as in rooftopMoM.basisFunctions.UniformRooftops: the x-directed rooftops first.

    properties (Access = protected)
        rooftops_ rooftopMoM.basisFunctions.UniformRooftops;
        spectra_ rooftopMoM.assembler.greensMatrix.Spectra;
        content_ rooftopMoM.basisFunctions.RooftopContent;
    end

    methods
        function obj = Applicator(rooftops, spectra, content)
            obj.rooftops_ = rooftops;
            obj.spectra_ = spectra;
            obj.content_ = content;
        end
    end

    methods
        function spectra = spectra(obj)
            spectra = obj.spectra_;
        end

        function bytes = memory(obj)
        % memory: Returns the bytes of the stored spectra.

            bytes = obj.spectra_.memory();
        end
    end

    methods (Abstract)
        % apply: Returns the product of the Green's matrix with an Nx1 vector of rooftop coefficients.
        output = apply(obj, input);

        % workingMemory: Returns the bytes of the arrays alive at the peak of one product, the output included.
        bytes = workingMemory(obj);
    end
end
