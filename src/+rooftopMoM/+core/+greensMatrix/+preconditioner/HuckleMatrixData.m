classdef HuckleMatrixData < handle
    properties (Access = private)
        context_;

        pxxf_;
        pxyf_;
        pyxf_;
        pyyf_;

        determinant_;

        rooftopContent_ rooftopMoM.basisFunctions.RooftopContent;
    end

    methods
        function obj = HuckleMatrixData(context)
            obj.context_ = context;
        end
    end

    methods
        function generate(obj, matrixData, pOption, centerScale)
            if nargin < 4
                centerScale = 1;
            end

            % Electric part is the largest and always present so we use it as an reference size.
            electric = matrixData.electric();
            [countX, countY] = size(matrixData.electric());

            magneticXX = rooftopMoM.core.greensMatrix.preconditioner.utils.padWithZeros(matrixData.magneticXX(), countX, countY);
            magneticYY = rooftopMoM.core.greensMatrix.preconditioner.utils.padWithZeros(matrixData.magneticYY(), countX, countY);

            switch pOption
                case 'full'
                    pX = countX; pY = countY;
                case 'half'
                    pX = floor(countX / 2); pY = floor(countY / 2);
                case 'sqrt'
                    pX = floor(sqrt(countX)); pY = floor(sqrt(countY));
                case 'one'
                    pX = 1; pY = 1;
                otherwise
                    error('Unknown pOption.');
            end

            magneticXXHuckle = rooftopMoM.core.greensMatrix.preconditioner.utils.generateHucklePreconditionerKernel(magneticXX, pX, pY, centerScale);
            magneticYYHuckle = rooftopMoM.core.greensMatrix.preconditioner.utils.generateHucklePreconditionerKernel(magneticYY, pX, pY, centerScale);
            electricHuckle = rooftopMoM.core.greensMatrix.preconditioner.utils.generateHucklePreconditionerKernel(electric, pX, pY, centerScale);

            rooftopContent = matrixData.rooftopContent();

            primitive = obj.context_.rooftops().primitive();
            [obj.pxxf_, obj.pxyf_, obj.pyxf_, obj.pyyf_] = rooftopMoM.core.greensMatrix.utils.generateSpectralGreensMatrix( ...
                magneticXXHuckle, magneticYYHuckle, electricHuckle, ...
                rooftopContent, ...
                primitive.sizeX(), primitive.sizeY() ...
            );

            obj.determinant_ = rooftopMoM.core.greensMatrix.preconditioner.utils.generateSpectralGreensMatrixDeterminant( ...
                obj.pxxf_, obj.pxyf_, obj.pyxf_, obj.pyyf_, ...
                rooftopContent() ...
            );

            obj.rooftopContent_ = rooftopContent;
        end
    end

    methods
        function pxxf = pxxf(obj)
            pxxf = obj.pxxf_;
        end

        function pxyf = pxyf(obj)
            pxyf = obj.pxyf_;
        end

        function pyxf = pyxf(obj)
            pyxf = obj.pyxf_;
        end

        function pyyf = pyyf(obj)
            pyyf = obj.pyyf_;
        end

        function determinant = determinant(obj)
            determinant = obj.determinant_;
        end

        function rooftopContent = rooftopContent(obj)
            rooftopContent = obj.rooftopContent_;
        end
    end
end
