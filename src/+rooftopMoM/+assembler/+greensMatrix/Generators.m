classdef Generators
% Generators: The three Toeplitz generators of the Green's matrix of a uniform grid.
%   magneticXX ((countX - 1) x countY) and magneticYY (countX x (countY - 1)) hold the vector potential part: the
%   interactions of the x- and y-directed rooftops with the reference rooftop of their direction. electric
%   (countX x countY) holds the scalar potential part: the interactions of the faces with the face at the origin.
%   The generator of a missing rooftop direction is empty.

    properties (Access = private)
        magneticXX_;
        magneticYY_;
        electric_;
        content_ rooftopMoM.basisFunctions.RooftopContent;
    end

    methods
        function obj = Generators(magneticXX, magneticYY, electric, content)
            obj.magneticXX_ = magneticXX;
            obj.magneticYY_ = magneticYY;
            obj.electric_ = electric;
            obj.content_ = content;
        end
    end

    methods
        function generators = modify(obj, modifier)
        % modify: Returns the generators with a function applied to each of them, e.g. @real or @imag for the real
        %   or the imaginary part of the Green's matrix.
        %
        %   Input:
        %       modifier: Function handle applied to every generator.
        %
        %   Output:
        %       generators: rooftopMoM.assembler.greensMatrix.Generators.

            generators = rooftopMoM.assembler.greensMatrix.Generators(modifier(obj.magneticXX_), ...
                modifier(obj.magneticYY_), modifier(obj.electric_), obj.content_);
        end

        function bytes = memory(obj)
            bytes = rooftopMoM.assembler.utils.generateByteCount(obj.magneticXX_, obj.magneticYY_, obj.electric_);
        end

        function magneticXX = magneticXX(obj)
            magneticXX = obj.magneticXX_;
        end

        function magneticYY = magneticYY(obj)
            magneticYY = obj.magneticYY_;
        end

        function electric = electric(obj)
            electric = obj.electric_;
        end

        function content = content(obj)
            content = obj.content_;
        end
    end
end
