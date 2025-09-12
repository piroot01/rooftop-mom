classdef MatrixData < handle
    properties (Access = private)
        data_ struct;

        rooftopContent_ rooftopMoM.basisFunctions.RooftopContent;

        paddingX_;
        paddingY_;
    end

    methods
        function obj = MatrixData(varargin)
            p = inputParser();

            addParameter(p, 'gxxFFT', []);
            addParameter(p, 'gxyFFT', []);
            addParameter(p, 'gyxFFT', []);
            addParameter(p, 'gyyFFT', []);

            addParameter(p, 'rooftopContent', rooftopMoM.basisFunctions.RooftopContent.NONE);
            addParameter(p, 'paddingX', []);
            addParameter(p, 'paddingY', []);

            parse(p, varargin{:});

            obj.data_ = struct( ...
                'gxxFFT', p.Results.gxxFFT, ...
                'gxyFFT', p.Results.gxyFFT, ...
                'gyxFFT', p.Results.gyxFFT, ...
                'gyyFFT', p.Results.gyyFFT ...
            );

            obj.rooftopContent_ = p.Results.rooftopContent;
            obj.paddingX_ = p.Results.paddingX;
            obj.paddingY_ = p.Results.paddingY;
        end
    end

    methods
        function data = data(obj)
            data = obj.data_;
        end

        function gxxFFT = gxxFFT(obj)
            gxxFFT = obj.data_.gxxFFT;
        end

        function gxyFFT = gxyFFT(obj)
            gxyFFT = obj.data_.gxyFFT;
        end

        function gyxFFT = gyxFFT(obj)
            gyxFFT = obj.data_.gyxFFT;
        end

        function gyyFFT = gyyFFT(obj)
            gyyFFT = obj.data_.gyyFFT;
        end

        function rooftopContent = rooftopContent(obj)
            rooftopContent = obj.rooftopContent_;
        end

        function paddingX = paddingX(obj)
            paddingX = obj.paddingX_;
        end

        function paddingY = paddingY(obj)
            paddingY = obj.paddingY_;
        end
    end
end
