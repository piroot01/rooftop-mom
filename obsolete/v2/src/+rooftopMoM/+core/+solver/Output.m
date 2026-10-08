classdef (Abstract) Output < handle
    properties (Access = private)
        name_ char = '';
        options_ struct = struct();

        % Information collected before the main loop.
        preRunInfo_ struct = struct();

        % Information collected after the main loop.
        postRunInfo_ struct = struct();
    end

    properties (Access = protected)
        iterationInfo_ cell;
    end

    methods
        function obj = Output(name, options)
            if nargin >= 1
                obj.name_ = name;
            end

            if nargin >= 2
                obj.options_ = options;
            end
        end
    end

    methods
        function preRunInfo = preRunInfo(obj)
            preRunInfo = obj.preRunInfo_;
        end

        function postRunInfo = postRunInfo(obj)
            postRunInfo = obj.postRunInfo_;
        end

        function iterationInfo = iterationInfo(obj)
            iterationInfo = obj.iterationInfo_;
        end

        function setPreRunInfo(obj, info)
            obj.preRunInfo_ = info;
        end

        function setPostRunInfo(obj, info)
            obj.postRunInfo_ = info;
        end

        function count = count(obj)
            count = numel(obj.iterationInfo_);
        end

        function value = selectIterationInfo(obj, fieldName, iterationNumber)
            if nargin < 3 || isempty(iterationNumber)
                value = cellfun(@(s) s.(fieldName), obj.iterationInfo_, 'UniformOutput', false);
            else
                if iterationNumber > obj.count()
                    error("The iteration number (%d) is out of range (1..%d).", iterationNumber, obj.count());
                else
                    value = obj.iterationInfo_{iterationNumber}.(fieldName);
                end
            end
        end
    end

    methods (Abstract)
        insertIterationInfo(obj, varargin);

        updateIterationInfo(obj, iteration, varargin);
    end
end
