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
            % Constructor: validate unit and store
            if ~isnumeric(value) || ~isscalar(value)
                error('Value must be a numeric scalar.');
            end

            if ~isfield(MemoryQuantity.UnitFactors, unit)
                error('Unsupported unit "%s".', unit);
            end

            obj.value_ = value;
            obj.unit_ = lower(unit);
        end
    end

    methods
        function bytes = bytes(obj)
            factor = obj.unitFactors_(obj.unit_);
            bytes = obj.value_ * factor;
        end
    end

    methods
        function memory = plus(lhs, rhs)
            assert(isa(lhs,'Memory') && isa(rhs,'Memory'), ...
                'Both operands must be Memory objects.');

            value = rhs.bytes + lhs.bytes;
            memory = Memory(value, 'B');
        end

        function memory = minus(lhs, rhs)
            assert(isa(lhs,'Memory') && isa(rhs,'Memory'), ...
                'Both operands must be Memory objects.');

            value = lhs.bytes - rhs.bytes;
            memory = Memory(value, 'B');
        end

        function memory = mtimes(lhs, rhs)
            assert(isa(lhs,'Memory') && isa(rhs,'Memory'), ...
                'Both operands must be Memory objects.');

            value = lhs.bytes * rhs.bytes;
            memory = Memory(value, 'B');
        end
    end
end
