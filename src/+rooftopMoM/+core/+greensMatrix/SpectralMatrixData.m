classdef SpectralMatrixData < handle
    properties (Access = private)
        context_;

        gxxf_;
        gxyf_;
        gyxf_;
        gyyf_;

        rooftopContent_ rooftopMoM.basisFunctions.RooftopContent;
    end

    methods
        function obj = SpectralMatrixData(context)
            obj.context_ = context;
        end
    end

    methods
        function generate(obj, matrixData)
            [magneticXXCirc, magneticYYCirc, electricCirc] = rooftopMoM.core.greensMatrix.utils.generateGreensMatrixCirculants(matrixData);

            primitive = obj.context_.rooftops().primitive();
            [obj.gxxf_, obj.gxyf_, obj.gyxf_, obj.gyyf_] = rooftopMoM.core.greensMatrix.utils.generateSpectralGreensMatrix( ...
                magneticXXCirc, magneticYYCirc, electricCirc, ...
                matrixData.rooftopContent(), ...
                primitive.sizeX(), primitive.sizeY() ...
            );

            obj.rooftopContent_ = matrixData.rooftopContent();
        end
    end

    methods
        function gxxf = gxxf(obj)
            gxxf = obj.gxxf_;
        end

        function gxyf = gxyf(obj)
            gxyf = obj.gxyf_;
        end

        function gyxf = gyxf(obj)
            gyxf = obj.gyxf_;
        end

        function gyyf = gyyf(obj)
            gyyf = obj.gyyf_;
        end

        function rooftopContent = rooftopContent(obj)
            rooftopContent = obj.rooftopContent_;
        end
    end
end
