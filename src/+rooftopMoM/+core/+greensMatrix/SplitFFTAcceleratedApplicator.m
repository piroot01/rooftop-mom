classdef SplitFFTAcceleratedApplicator < rooftopMoM.core.greensMatrix.GreensMatrixApplicator
    properties (Access = private)
        gxxf_;
        gxyf_;
        gyxf_;
        gyyf_;

        determinant_;

        % phase-shift arrays
        pr_;
        pc_;
        prc_;
        pcc_;
    end

    methods
        function obj = SplitFFTAcceleratedApplicator(context, matrixData)
            obj@rooftopMoM.core.greensMatrix.GreensMatrixApplicator(context, matrixData);
        end
    end

    methods
        function generate(obj)
            primitive = obj.context_.rooftops().primitive();
            [gridSizeX, gridSizeY] = primitive.gridSize();

            [magneticXXCirc, magneticYYCirc, electricCirc] = rooftopMoM.core.greensMatrix.utils.generateGreensMatrixCirculants(obj.matrixData_, gridSizeX + 1, gridSizeY + 1);

            [gxxf, gxyf, gyxf, gyyf, obj.determinant_] = rooftopMoM.core.greensMatrix.utils.generateSpectralGreensMatrix( ...
                magneticXXCirc, magneticYYCirc, electricCirc, ...
                obj.matrixData_.rooftopContent(), ...
                primitive.sizeX(), ...
                primitive.sizeY() ...
            );

            % Subsample each of the matrices.
            obj.gxxf_ = rooftopMoM.core.greensMatrix.utils.subsampleMatrix(gxxf);
            clear gxxf;

            obj.gxyf_ = rooftopMoM.core.greensMatrix.utils.subsampleMatrix(gxyf);
            clear gxyf;

            obj.gyxf_ = rooftopMoM.core.greensMatrix.utils.subsampleMatrix(gyxf);
            clear gyxf;

            obj.gyyf_ = rooftopMoM.core.greensMatrix.utils.subsampleMatrix(gyyf);
            clear gyyf;

            M = double(gridSizeX);
            N = double(gridSizeY);

            % Generate the phase-shift arrays.
            kr = (0:(M - 1)).';
            kc = (0:(N - 1));

            obj.pr_ = exp(-1i * pi * kr / M);
            obj.pc_ = exp(-1i * pi * kc / N);
            obj.prc_ = conj(obj.pr_);
            obj.pcc_ = conj(obj.pc_);
        end

        function output = apply(obj, input)
            primitive = obj.context_.rooftops().primitive();
            countX = primitive.countX();
            countY = primitive.countY();

            % all blocks of the greens matrix have same size
            sizeX = obj.context_.rooftops().countDir('x');
            sizeY = obj.context_.rooftops().countDir('y');

            switch obj.matrixData_.rooftopContent()
                case rooftopMoM.basisFunctions.RooftopContent.XONLY
                    inputXReshaped = zeros(countX, countY);
                    inputXReshaped(1:(end - 1)) = input;

                    inputXReshapedSpectr_e = fft(inputXReshaped, [], 1);
                    inputXReshapedSpectr_o = fft(inputXReshaped .* obj.pr_, [], 1);
                    clear inputXReshaped;

                    % outputX = gxxf .* inputX
                    outputXSpectr_ee = obj.gxxf_.ee .* inputXReshapedSpectr_e;
                    outputXSpectr_eo = obj.gxxf_.eo .* inputXReshapedSpectr_e;
                    clear inputXReshapedSpectr_e;
                    outputXSpectr_oe = obj.gxxf_.oe .* inputXReshapedSpectr_o;
                    outputXSpectr_oo = obj.gxxf_.oo .* inputXReshapedSpectr_o;
                    clear inputXReshapedSpectr_o;

                    outputXSpectr_e = outputXSpectr_ee + outputXSpectr_eo;
                    clear outputXSpectr_ee outputXSpectr_eo;

                    outputXSpectr_o = outputXSpectr_oe + outputXSpectr_oo;
                    clear outputXSpectr_oe outputXSpectr_oo;

                    outputX = 0.25 * (ifft(outputXSpectr_e, [], 1) + ifft(outputXSpectr_o, [], 1) .* obj.prc_);
                    clear outputXSpectr_e outputXSpectr_o;

                    output = outputX(1:sizeX);
                case rooftopMoM.basisFunctions.RooftopContent.YONLY
                    input = input.';
                    inputYReshaped = zeros(countX, countY);
                    inputYReshaped(1:(end - 1)) = input;

                    inputYReshapedSpectr_e = fft(inputYReshaped, [], 2);
                    inputYReshapedSpectr_o = fft(inputYReshaped .* obj.pc_, [], 2);

                    outputYSpectr_ee = obj.gyyf_.ee .* inputYReshapedSpectr_e;
                    outputYSpectr_oe = obj.gyyf_.oe .* inputYReshapedSpectr_e;
                    clear inputYReshapedSpectr_e;
                    outputYSpectr_eo = obj.gyyf_.eo .* inputYReshapedSpectr_o;
                    outputYSpectr_oo = obj.gyyf_.oo .* inputYReshapedSpectr_o;
                    clear inputYReshapedSpectr_o;

                    outputYSpectr_e = outputYSpectr_ee + outputYSpectr_oe;
                    clear outputYSpectr_ee outputYSpectr_oe;

                    outputYSpectr_o = outputYSpectr_eo + outputYSpectr_oo;
                    clear outputYSpectr_eo outputYSpectr_oo;

                    outputY = 0.25 * (ifft(outputYSpectr_e, [], 2) + ifft(outputYSpectr_o, [], 2) .* obj.pcc_);
                    clear outputYSpectr_e outputYSpectr_o;

                    output = outputY(1:sizeY).';
                case rooftopMoM.basisFunctions.RooftopContent.BOTH
                    % First we process the x direction of the input
                    % We extend the input in x direction such that the size matches the subsampled greens matrix
                    inputXReshaped = zeros(countX, countY);
                    inputXReshaped(1:(end - 1), 1:end) = reshape(input(1:sizeX), countX - 1, countY);

                    inputXReshapedSpectr_e = fft(inputXReshaped, [], 1);

                    inputXReshapedSpectr_ee = fft(inputXReshapedSpectr_e, [], 2);
                    inputXReshapedSpectr_eo = fft(inputXReshapedSpectr_e .* obj.pc_, [], 2);
                    clear inputXReshapedSpectr_e;

                    inputXReshapedSpectr_o = fft(inputXReshaped .* obj.pr_, [], 1);
                    clear inputXReshaped;

                    inputXReshapedSpectr_oe = fft(inputXReshapedSpectr_o, [], 2);
                    inputXReshapedSpectr_oo = fft(inputXReshapedSpectr_o .* obj.pc_, [], 2);
                    clear inputXReshapedSpectr_o;

                    % outputX = gxxf .* inputX
                    outputXSpectr_ee = obj.gxxf_.ee .* inputXReshapedSpectr_ee;
                    outputXSpectr_eo = obj.gxxf_.eo .* inputXReshapedSpectr_eo;
                    outputXSpectr_oe = obj.gxxf_.oe .* inputXReshapedSpectr_oe;
                    outputXSpectr_oo = obj.gxxf_.oo .* inputXReshapedSpectr_oo;

                    % outputY = gyxf .* inputX
                    outputYSpectr_ee = obj.gyxf_.ee .* inputXReshapedSpectr_ee;
                    clear inputXReshapedSpectr_ee;
                    outputYSpectr_eo = obj.gyxf_.eo .* inputXReshapedSpectr_eo;
                    clear inputXReshapedSpectr_eo;
                    outputYSpectr_oe = obj.gyxf_.oe .* inputXReshapedSpectr_oe;
                    clear inputXReshapedSpectr_oe;
                    outputYSpectr_oo = obj.gyxf_.oo .* inputXReshapedSpectr_oo;
                    clear inputXReshapedSpectr_oo;

                    % Secondly we process the y direction of the input
                    % The same logic as in x direction case
                    inputYReshaped = zeros(countX, countY);
                    inputYReshaped(1:end, 1:(end - 1)) = reshape(input((sizeX + 1):end), countX, countY - 1);

                    inputYReshapedSpectr_e = fft(inputYReshaped, [], 1);

                    inputYReshapedSpectr_ee = fft(inputYReshapedSpectr_e, [], 2);
                    inputYReshapedSpectr_eo = fft(inputYReshapedSpectr_e .* obj.pc_, [], 2);
                    clear inputYReshapedSpectr_e;

                    inputYReshapedSpectr_o = fft(inputYReshaped .* obj.pr_, [], 1);
                    clear inputYReshaped;

                    inputYReshapedSpectr_oe = fft(inputYReshapedSpectr_o, [], 2);
                    inputYReshapedSpectr_oo = fft(inputYReshapedSpectr_o .* obj.pc_, [], 2);
                    clear inputYReshapedSpectr_o;

                    % outputX = outputX + gxyf .* inputX;
                    outputXSpectr_ee = outputXSpectr_ee + obj.gxyf_.ee .* inputYReshapedSpectr_ee;
                    outputXSpectr_eo = outputXSpectr_eo + obj.gxyf_.eo .* inputYReshapedSpectr_eo;
                    outputXSpectr_oe = outputXSpectr_oe + obj.gxyf_.oe .* inputYReshapedSpectr_oe;
                    outputXSpectr_oo = outputXSpectr_oo + obj.gxyf_.oo .* inputYReshapedSpectr_oo;

                    % outputY = gyyf .* inputY
                    outputYSpectr_ee = outputYSpectr_ee + obj.gyyf_.ee .* inputYReshapedSpectr_ee;
                    clear inputYReshapedSpectr_ee;
                    outputYSpectr_eo = outputYSpectr_eo + obj.gyyf_.eo .* inputYReshapedSpectr_eo;
                    clear inputYReshapedSpectr_eo;
                    outputYSpectr_oe = outputYSpectr_oe + obj.gyyf_.oe .* inputYReshapedSpectr_oe;
                    clear inputYReshapedSpectr_oe;
                    outputYSpectr_oo = outputYSpectr_oo + obj.gyyf_.oo .* inputYReshapedSpectr_oo;
                    clear inputYReshapedSpectr_oo;

                    outputXSpectr_e = ifft(outputXSpectr_ee, [], 2) + ifft(outputXSpectr_eo, [], 2) .* obj.pcc_;
                    clear outputXSpectr_ee outputXSpectr_eo;

                    outputXSpectr_o = ifft(outputXSpectr_oe, [], 2) + ifft(outputXSpectr_oo, [], 2) .* obj.pcc_;
                    clear outputXSpectr_oe outputXSpectr_oo;

                    outputX = 0.25 * (ifft(outputXSpectr_e, [], 1) + ifft(outputXSpectr_o, [], 1) .* obj.prc_);
                    clear outputXSpectr_e outputXSpectr_o;

                    outputYSpectr_e = ifft(outputYSpectr_ee, [], 2) + ifft(outputYSpectr_eo, [], 2) .* obj.pcc_;
                    clear outputYSpectr_ee outputYSpectr_eo;

                    outputYSpectr_o = ifft(outputYSpectr_oe, [], 2) + ifft(outputYSpectr_oo, [], 2) .* obj.pcc_;
                    clear outputYSpectr_oe outputYSpectr_oo;

                    outputY = 0.25 * (ifft(outputYSpectr_e, [], 1) + ifft(outputYSpectr_o, [], 1) .* obj.prc_);
                    clear outputYSpectr_e outputYSpectr_o;

                    % Reshape back to column vectors
                    outputX = outputX(1:(end - 1), 1:end);
                    outputY = outputY(1:end, 1:(end - 1));

                    output = [outputX(:); outputY(:)];
                otherwise
                    error('Unknown rooftop content: %s', rooftopContent);
            end
        end

        function varargout = inputSize(obj)
            if nargout <= 1
                varargout{1} = obj.context_.rooftops().count();
            elseif nargout == 2
                varargout{1} = obj.context_.rooftops().count();
                varargout{2} = 1;
            end
        end

        function memory = memory(obj)
            % TODO: implement
        end
    end
end
