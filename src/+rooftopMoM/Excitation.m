classdef Excitation
% Excitation: Right-hand side of a Problem, an incident field tested with every rooftop. It is created by one of the
%   sources:
%       excitation = rooftopMoM.Excitation.planeWave(problem, 'direction', [0; 0; -1], 'polarization', [1; 0; 0]);
%       excitation = rooftopMoM.Excitation.deltaGap(problem, 'voltage', 1);

    properties (Access = private)
        problem_ rooftopMoM.Problem;

        % 'planeWave' or 'deltaGap'.
        type_ char;

        vector_ double;

        % Struct of the source: direction, polarization and intensity of a plane wave, port and voltage of a delta gap.
        parameters_ struct;
    end

    methods (Access = private)
        function obj = Excitation(problem, type, vector, parameters)
            obj.problem_ = problem;
            obj.type_ = type;
            obj.vector_ = vector;
            obj.parameters_ = parameters;
        end
    end

    methods (Static)
        function excitation = planeWave(problem, options)
        % planeWave: Returns the excitation of an incident plane wave,
        %   E(r) = intensity * polarization * exp(-1i * wavenumber * direction.' * r).
        %
        %   Input:
        %       problem: rooftopMoM.Problem.
        %
        %   Optional name-value inputs:
        %       'direction': 3x1 direction of propagation (default [0; 0; -1], normal incidence from above).
        %       'polarization': 3x1 direction of the electric field, perpendicular to the direction (default
        %           [1; 0; 0]).
        %       'intensity': Amplitude of the electric field in V/m (default 1).
        %
        %   Output:
        %       excitation: rooftopMoM.Excitation.

            arguments
                problem (1, 1) rooftopMoM.Problem
                options.direction (3, 1) {mustBeReal, mustBeFinite} = [0; 0; -1]
                options.polarization (3, 1) {mustBeReal, mustBeFinite} = [1; 0; 0]
                options.intensity (1, 1) {mustBeReal, mustBeFinite, mustBePositive} = 1
            end

            direction = normalize_(double(options.direction), 'direction');
            polarization = normalize_(double(options.polarization), 'polarization');

            % Normalizing and summing round the dot product of perpendicular unit vectors by at most 5 * eps; twice
            % that also covers inputs that are rounded themselves.
            dotProduct = direction.' * polarization;

            if abs(dotProduct) > 10 * eps
                error('rooftopMoM:excitation:notPerpendicular', ['The polarization must be perpendicular to the ' ...
                    'direction of propagation, their dot product is %.1e.'], dotProduct);
            end

            parameters = struct('direction', direction, 'polarization', polarization, ...
                'intensity', double(options.intensity));
            vector = rooftopMoM.operator.generatePlaneWave(problem.rooftops(), problem.wavenumber(), direction, ...
                polarization, parameters.intensity);
            excitation = rooftopMoM.Excitation(problem, 'planeWave', vector, parameters);

            rooftopMoM.logging.Logger.instance().debug(['Plane wave in the direction [%g %g %g], polarized along ' ...
                '[%g %g %g], %g V/m.'], direction, polarization, parameters.intensity);
        end

        function excitation = deltaGap(problem, options)
        % deltaGap: Returns the excitation of a delta-gap voltage source across the edges of a port.
        %
        %   Input:
        %       problem: rooftopMoM.Problem.
        %
        %   Optional name-value inputs:
        %       'port': Indices of the rooftops across the gap (default the rooftops of the given direction on the
        %           middle line of edges, see rooftopMoM.operator.utils.generateMiddlePort).
        %       'direction': 'x' or 'y', the direction of the default port (default 'x').
        %       'voltage': Voltage of the source in V (default 1).
        %
        %   Output:
        %       excitation: rooftopMoM.Excitation.

            arguments
                problem (1, 1) rooftopMoM.Problem
                options.port (1, :) {mustBeInteger, mustBePositive} = []
                options.direction (1, 1) char {mustBeMember(options.direction, {'x', 'y'})} = 'x'
                options.voltage (1, 1) {mustBeNumeric, mustBeFinite} = 1
            end

            rooftops = problem.rooftops();
            port = double(options.port);

            if isempty(port)
                if rooftops.countDir(options.direction) == 0
                    error('rooftopMoM:excitation:noRooftops', 'The mesh has no %s-directed rooftops.', ...
                        options.direction);
                end

                port = rooftopMoM.operator.utils.generateMiddlePort(rooftops, options.direction);
            elseif any(port > rooftops.count()) || numel(unique(port)) < numel(port)
                error('rooftopMoM:excitation:invalidPort', 'The port must hold distinct rooftop indices up to %d.', ...
                    rooftops.count());
            end

            parameters = struct('port', port, 'voltage', double(options.voltage));
            vector = rooftopMoM.operator.generateDeltaGap(rooftops, port, parameters.voltage);
            excitation = rooftopMoM.Excitation(problem, 'deltaGap', vector, parameters);

            rooftopMoM.logging.Logger.instance().debug('Delta gap of %g V across %d rooftops.', ...
                abs(parameters.voltage), numel(port));
        end
    end

    methods
        function problem = problem(obj)
            problem = obj.problem_;
        end

        function type = type(obj)
        % type: Returns the source, 'planeWave' or 'deltaGap'.

            type = obj.type_;
        end

        function vector = vector(obj)
        % vector: Returns the Nx1 excitation vector.

            vector = obj.vector_;
        end

        function parameters = parameters(obj)
        % parameters: Returns the struct of the source: the fields direction, polarization and intensity of a plane
        %   wave, the fields port and voltage of a delta gap.

            parameters = obj.parameters_;
        end
    end
end

function unit = normalize_(vector, name)
    if ~any(vector)
        error('rooftopMoM:excitation:zeroVector', 'The %s must not be zero.', name);
    end

    unit = vector / norm(vector);
end
