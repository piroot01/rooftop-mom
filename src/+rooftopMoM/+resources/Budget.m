classdef Budget < handle
% Budget: Shared memory budget of the objects of rooftopMoM. It keeps the ceiling and the chunk size of the assembly
%   given by the user and the bytes held by the live objects. New memory is admitted only if it fits below the ceiling
%   and into the free memory of the machine.

    properties (Constant, Access = private)
        % Chunks of the assembly of about the size of the last-level cache are the fastest; larger ones are slower.
        defaultChunkSize_ = 2^24;
    end

    properties (Access = private)
        % Bytes the objects may hold together, Inf for no ceiling.
        ceiling_ double = Inf;

        % Bytes of one chunk of the assembly given by the user, NaN for the default.
        chunkSize_ double = NaN;

        % Bytes held by the live objects.
        held_ double = 0;
    end

    methods (Access = private)
        function obj = Budget()
        end
    end

    methods (Static)
        function budget = instance()
        % instance: Returns the budget shared by all objects.

            persistent shared;

            if isempty(shared) || ~isvalid(shared)
                shared = rooftopMoM.resources.Budget();
            end

            budget = shared;
        end
    end

    methods
        function configure(obj, ceiling, chunkSize)
        % configure: Sets the ceiling in bytes, Inf for none, and the chunk size of the assembly in bytes, NaN for the
        %   default.

            obj.ceiling_ = ceiling;
            obj.chunkSize_ = chunkSize;
        end

        function ceiling = ceiling(obj)
            ceiling = obj.ceiling_;
        end

        function held = held(obj)
            held = obj.held_;
        end

        function hold(obj, bytes)
        % hold: Adds the bytes of a new object to the held memory.

            obj.held_ = obj.held_ + bytes;
        end

        function release(obj, bytes)
        % release: Removes the bytes of a deleted object from the held memory.

            obj.held_ = obj.held_ - bytes;
        end

        function bytes = chooseChunk(obj, otherBytes)
        % chooseChunk: Returns the bytes of one chunk of the assembly: the chunk size given by the user as it is, or
        %   else the default, at most what the ceiling leaves next to otherBytes.

            bytes = obj.chunkSize_;

            if isnan(bytes)
                bytes = min(obj.defaultChunkSize_, max(0, obj.ceiling_ - obj.held_ - otherBytes));
            end
        end

        function require(obj, bytes, description)
        % require: Raises an error if the bytes do not fit below the ceiling or into the free memory of the machine.
        %
        %   Inputs:
        %       bytes: Memory the new data needs.
        %       description: Character vector naming the data, the subject of the message.

            available = rooftopMoM.resources.readAvailableMemory();
            toText = @rooftopMoM.resources.utils.formatBytes;

            rooftopMoM.logging.Logger.instance().debug('%s needs %s; held %s, ceiling %s, free %s.', description, ...
                toText(bytes), toText(obj.held_), toText(obj.ceiling_), toText(available));

            if obj.held_ + bytes > obj.ceiling_
                error('rooftopMoM:resources:insufficientMemory', ...
                    '%s needs %s, but only %s of the memory ceiling of %s are left.', description, toText(bytes), ...
                    toText(max(0, obj.ceiling_ - obj.held_)), toText(obj.ceiling_));
            end

            % An unknown free memory leaves the ceiling as the only limit.
            if bytes > available
                error('rooftopMoM:resources:insufficientMemory', '%s needs %s, but only %s are free on this machine.', ...
                    description, toText(bytes), toText(available));
            end
        end
    end
end
