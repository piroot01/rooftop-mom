classdef DeltaGapFeederExcitation < rooftopMoM.operators.OperatorInterface
    methods (Static)
        function V = generate(basisFunctions, mode, voltage)
        % generate: Build a voltage vector from a binary mode mask.
        %
        %   Inputs:
        %       basisFunctions: Object implementing BasisFunctionsInterface.
        %           Must provide count() for number of ports (N).
        %       mode: Nx1 logical array indicating active ports.
        %       voltage: Scalar or Mx1 numeric vector, where M = sum(mode).
        %           If scalar, same voltage applied to all active ports.
        %           If vector, each element maps to one active port.
        %
        %   Output:
        %       V: Nx1 numeric vector, with V(i)=voltage(k) for the k-th true element in mode, otherwise 0.

            V = zeros(basisFunctions.count(), 1);

            % Assign the specified voltage to each active port index
            V(mode) = voltage;
        end
    end
end
