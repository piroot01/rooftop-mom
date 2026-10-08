classdef (Abstract) Base < unit.Base
% Base: Common fixture of the solver tests, small dense systems of a known condition.

    methods (TestMethodSetup)
        function silenceLogger(~)
            % Some tests do not converge on purpose; their warnings would clutter the test output.
            rooftopMoM.logging.Logger.instance().configure(rooftopMoM.logging.LogLevel.OFF, {});
        end
    end

    methods (Access = protected)
        function [matrix, rightHandSide] = createSymmetricSystem(~, count)
        % createSymmetricSystem: Returns a complex symmetric system, 4 * I plus a symmetric random part of norm about 2,
        %   so its singular values lie between about 2 and 6.

            rng(1);
            random = complex(randn(count), randn(count)) / sqrt(count);
            matrix = 4 * eye(count) + (random + random.') / 2;
            rightHandSide = complex(randn(count, 1), randn(count, 1));
        end

        function [matrix, rightHandSide] = createNonsymmetricSystem(~, count)
        % createNonsymmetricSystem: Returns a general system, 4 * I plus a random part with its eigenvalues within
        %   about 1.4 of zero.

            rng(2);
            matrix = 4 * eye(count) + complex(randn(count), randn(count)) / sqrt(count);
            rightHandSide = complex(randn(count, 1), randn(count, 1));
        end
    end
end
