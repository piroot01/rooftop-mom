classdef MatrixData < handle
    properties (Access = private)
        magneticXX_;
        magneticYY_;
        electric_;

        rooftopContent_ rooftopMoM.basisFunctions.RooftopContent;
    end

    methods
        function obj = MatrixData(varargin)
            p = inputParser();

            addParameter(p, 'magneticXX', []);
            addParameter(p, 'magneticYY', []);
            addParameter(p, 'electric', []);

            addParameter(p, 'rooftopContent', rooftopMoM.basisFunctions.RooftopContent.NONE);

            parse(p, varargin{:});

            obj.magneticXX_ = p.Results.magneticXX;
            obj.magneticYY_ = p.Results.magneticYY;
            obj.electric_ = p.Results.electric;

            obj.rooftopContent_ = p.Results.rooftopContent;
        end
    end

    methods
        function magneticXX = magneticXX(obj)
            magneticXX = obj.magneticXX_;
        end

        function magneticYY = magneticYY(obj)
            magneticYY = obj.magneticYY_;
        end

        function electric = electric(obj)
            electric = obj.electric_;
        end

        function rooftopContent = rooftopContent(obj)
            rooftopContent = obj.rooftopContent_;
        end
    end
end
