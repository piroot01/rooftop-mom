classdef InputImpedance < rooftopMoM.operators.OperatorInterface
    methods (Static)
        function Zin = generate(impedanceMatrix, currentVector, mode)
        % generate: Compute the input impedance seen at the excited port(s).
        %
        %   Inputs:
        %       impedanceMatrix: N×N impedance matrix Z.
        %       currentVector: N×1 vector I of port currents.
        %       mode: Nx1 logical array indicating active ports.
        %
        %   Output:
        %       Zin: scalar (or vector) of input impedance(s),
        %           Z_in = I^H Z I ./ |I_f|^2,
        %           where I_f is the current at the excited port(s).

            % numerator: I^H * Z * I
            numerator = currentVector' * impedanceMatrix * currentVector;

            % denominator: squared magnitude of the excitation current(s).
            denom = abs(currentVector(mode)).^2;

            Zin = numerator ./ denom;
        end
    end
end
