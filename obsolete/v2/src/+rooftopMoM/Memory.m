classdef Memory < handle
    properties (Access = private)
        value_ double;
        unit_ char;
    end

    properties (Constant, Access = private)
        % Map of supported units to their factor in bytes
        unitFactors_ = struct( ...
            'b', 1, ...
            'kb', 1e3, ...
            'mb', 1e6, ...
            'gb', 1e9, ...
            'tb', 1e12, ...
            'kib', 2^10, ...
            'mib', 2^20, ...
            'gib', 2^30, ...
            'tib', 2^40 ...
        );
    end

    methods
        function obj = Memory(value, unit)
            unit = lower(unit);

            % Constructor: validate unit and store
            if ~isnumeric(value) || ~isscalar(value)
                error('Value must be a numeric scalar.');
            end

            if ~isfield(rooftopMoM.Memory.unitFactors_, unit)
                error('Unsupported unit "%s".', unit);
            end

            obj.value_ = value;
            obj.unit_ = unit;
        end
    end

    methods
        function bytes = bytes(obj)
            factor = obj.unitFactors_.(obj.unit_);
            bytes = obj.value_ * factor;
        end
    end

    methods
        function memory = plus(lhs, rhs)
            assert(isa(lhs,'rooftopMoM.Memory') && isa(rhs,'rooftopMoM.Memory'), ...
                'Both operands must be rooftopMoM.Memory objects.');

            value = rhs.bytes + lhs.bytes;
            memory = rooftopMoM.Memory(value, 'B');
        end

        function memory = minus(lhs, rhs)
            assert(isa(lhs,'rooftopMoM.Memory') && isa(rhs,'rooftopMoM.Memory'), ...
                'Both operands must be rooftopMoM.Memory objects.');

            value = lhs.bytes - rhs.bytes;
            memory = rooftopMoM.Memory(value, 'B');
        end

        function memory = mtimes(lhs, rhs)
            % Support both scalar*Memory and Memory*scalar
            if isa(lhs,'rooftopMoM.Memory') && isnumeric(rhs) && isscalar(rhs)
                bytes = lhs.bytes() * rhs;
            elseif isa(rhs,'rooftopMoM.Memory') && isnumeric(lhs) && isscalar(lhs)
                bytes = rhs.bytes() * lhs;
            else
                error('One operand must be a scalar and the other a rooftopMoM.Memory.');
            end

            % Always return a Memory object in bytes
            memory = rooftopMoM.Memory(bytes, 'B');
        end
    end
end
