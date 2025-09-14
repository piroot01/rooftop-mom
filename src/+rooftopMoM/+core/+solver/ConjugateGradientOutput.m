classdef ConjugateGradientOutput < rooftopMoM.core.solver.Output
    methods
        function obj = ConjugateGradientOutput(name, options)
            if nargin < 1
                name = 'ConjugateGradientOutput';
            end

            if nargin < 2
                options = struct();
            end

            obj@rooftopMoM.core.solver.Output(name, options);
        end
    end

    methods
        function insertIterationInfo(obj, varargin)
            p = inputParser;

            addParameter(p, 'residueNorm', []);

            parse(p, varargin{:});

            iterationInfo = struct( ...
                'residueNorm', p.Results.residueNorm ...
            );

            obj.iterationInfo_{end + 1} = iterationInfo;
        end

        function updateIterationInfo(obj, iterationNumber, varargin)
        % updateIterationInfo: Update selected fields of an existing iteration entry.
        %   This function accepts the same parameter names as insertInterationInfo.

            total = obj.count();
            if iterationNumber < 1 || iterationNumber > total
                error("The iteration number (%d) is out of range (1..%d).", iterationNumber, total);
            end

            % parse provided name/value pairs (same defaults as insertInterationInfo)
            p = inputParser;
            addParameter(p, 'residueNorm', []);

            parse(p, varargin{:});

            % Determine which parameters were actually provided (not defaults)
            providedParams = setdiff(fieldnames(p.Results), p.UsingDefaults);

            % nothing to update -> do nothing
            if isempty(providedParams)
                return;
            end

            entry = obj.iterationInfo_{iterationNumber};

            for k = 1:numel(providedParams)
                fn = providedParams{k};
                entry.(fn) = p.Results.(fn);
            end
            
            obj.iterationInfo_{iterationNumber} = entry;
        end
    end
end
