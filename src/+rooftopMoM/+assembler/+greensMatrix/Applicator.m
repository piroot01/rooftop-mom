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

        function bytes = workingMemory(obj)
        % workingMemory: Returns the bytes of the arrays alive at the peak of one product, the output included.

            [sizeX, sizeY] = obj.spectra_(1).gridSize();
            bytes = obj.estimateWorkingMemory(obj.rooftops_, sizeX, sizeY, obj.content_);
        end
    end

    methods (Abstract)
        % apply: Returns the product of the Green's matrix with an Nx1 vector of rooftop coefficients.
        output = apply(obj, input);
    end

    methods (Abstract, Static)
        % estimateWorkingMemory: Returns the bytes of the arrays alive at the peak of one product, the output included,
        %   for the spectra on a sizeX x sizeY grid, before any spectra exist.
        bytes = estimateWorkingMemory(rooftops, sizeX, sizeY, content);
    end
end
