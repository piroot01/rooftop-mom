classdef Result
% Result: Outcome of rooftopMoM.Solver.solve: the solution and its convergence record. For an Excitation the
%   solution holds the rooftop coefficients of the current, and the result gives the input impedance of a delta gap
%   or the radar cross section of a plane wave.

    properties (Access = private)
        output_ rooftopMoM.solver.Output;

        % rooftopMoM.Excitation of the right-hand side, empty for a plain vector.
        excitation_;
    end

    methods
        function obj = Result(output, excitation)
            obj.output_ = output;
            obj.excitation_ = excitation;
        end
    end

    methods
        function solution = solution(obj)
            solution = obj.output_.solution();
        end

        function converged = converged(obj)
            converged = obj.output_.converged();
        end

        function iterationCount = iterationCount(obj)
            iterationCount = obj.output_.iterationCount();
        end

        function residualNorms = residualNorms(obj)
        % residualNorms: Returns the relative residual norms monitored by the method, from the initial guess to the
        %   last iteration.

            residualNorms = obj.output_.residualNorms();
        end

        function residualNorm = residualNorm(obj)
        % residualNorm: Returns the relative norm of the true residual, norm(b - A * x) / norm(b).

            residualNorm = obj.output_.residualNorm();
        end

        function time = time(obj)
        % time: Returns the time of the solve in seconds.

            time = obj.output_.time();
        end

        function excitation = excitation(obj)
            excitation = obj.excitation_;
        end

        function impedance = inputImpedance(obj)
        % inputImpedance: Returns the input impedance of the delta-gap source in ohms.

            parameters = obj.parameters_('deltaGap');
            impedance = rooftopMoM.operator.computeInputImpedance(obj.excitation_.problem().rooftops(), ...
                obj.output_.solution(), parameters.port, parameters.voltage);
        end

        function crossSection = radarCrossSection(obj)
        % radarCrossSection: Returns the monostatic radar cross section of the plane wave in m^2, for its
        %   polarization.

            parameters = obj.parameters_('planeWave');
            crossSection = rooftopMoM.operator.computeRadarCrossSection(obj.excitation_.problem().wavenumber(), ...
                obj.excitation_.vector(), obj.output_.solution(), parameters.intensity);
        end
    end

    methods (Access = private)
        function parameters = parameters_(obj, type)
        % parameters_: Returns the parameters of the source, which must be of the given type.

            if isempty(obj.excitation_) || ~strcmp(obj.excitation_.type(), type)
                error('rooftopMoM:result:wrongExcitation', 'This quantity needs a %s excitation.', type);
            end

            parameters = obj.excitation_.parameters();
        end
    end
end
