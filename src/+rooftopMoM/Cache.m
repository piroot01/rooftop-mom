classdef Cache < handle
% Cache: Singleton cache manager for rooftopMoM data.
%
%   NOTE:
%       - Singleton instance is private and managed by getInstance_().
%       - By default global caching is disabled and no options are enabled.

    properties (Access = private)
        % Nested struct is used for holding cached data.
        store_;

        % containers.Map string->logical for option switches (key: dot-path)
        options_;

        cacheEnabled_ logical;
    end

    methods (Access = private)
        function obj = Cache()
            obj.store_ = struct();
            obj.options_ = containers.Map('KeyType','char','ValueType','any');
            obj.cacheEnabled_ = false;
        end
    end

    methods (Static, Access = private)
        function instance = getInstance_()
            persistent INST

            if isempty(INST) || ~isvalid(INST)
                INST = rooftopMoM.Cache();
            end

            instance = INST;
        end
    end

     methods (Access = private)
        % Global enable/disable
        function setStatus_(obj, flag)
            obj.cacheEnabled_ = logical(flag);
        end

        function flag = getStatus_(obj)
            flag = logical(obj.cacheEnabled_);
        end

        % Options management
        function setOption_(obj, path, flag)
            obj.options_(char(path)) = logical(flag);
        end

        function setOptionFromKey_(obj, key)
            key = char(key);

            if endsWith(k, '_enable')
                base = key(1:(end - 7));
                flag = true;
            elseif endsWith(k, '_disable')
                base = key(1:(end - 8));
                flag = false;
            else
                error('Invalid key format!');
            end

            path = strrep(base, '_', '.');
            obj.setOption_(path, flag);
        end

        function tf = getOptionStatus_(obj, path)
            path = char(path);

            if obj.options_.isKey(path)
                tf = logical(obj.options_(path));
            else
                tf = false;
            end
        end

        function keys = listOptions_(obj)
            keys = obj.options_.keys();
        end

        % Storage API
        function autoSet_(obj, path, data)
        % Check global, option and set only if both are enabled

            if obj.getStatus_() && obj.getOptionStatus_(path)
                obj.set_(path, data);
            end
        end

        function set_(obj, path, data)
        % Store data at nested dot-path regardless of option (use autoSet_ to gate)

            parts = strsplit(char(path), '.');
            obj.store_ = rooftopMoM.Cache.setNested_(obj.store_, parts, data);
        end

        function data = get_(obj, path)
            parts = strsplit(char(path), '.');
            [flag, data] = rooftopMoM.Cache.getNested_(obj.store_, parts);

            if ~flag
                data = [];
            end
        end

        function cleared = clear_(obj, path)
        % Clear subtree at path; if path omitted clear whole store.

            if nargin < 2 || isempty(path)
                obj.store_ = struct();
                cleared = true;
                return;
            end

            parts = strsplit(char(path), '.');
            obj.store_ = rooftopMoM.Cache.clearNested_(obj.store_, parts);
            cleared = true;
        end

        function snapshot = snapshot_(obj)
            snapshot = obj.store_;
        end

        function destroy_(obj)
        % Reset data and options and global flag.

            obj.store_ = struct();
            obj.options_ = containers.Map('KeyType','char','ValueType','any');
            obj.cacheEnabled_ = false;
            % NOTE: getInstance_() persistent INST will remain; instance is now reset.
        end
    end

    % Private static helpers for nested struct manipulation
    methods (Static, Access = private)
        function store = setNested_(store, parts, data)
            if isempty(parts)
                store = data;
                return;
            end

            part = parts{1};

            if ~isfield(store, part) || ~isstruct(store.(part))
                store.(part) = struct();
            end

            store.(part) = rooftopMoM.Cache.setNested_(store.(part), parts(2:end), data);
        end

        function [found, data] = getNested_(store, parts)
            if isempty(parts)
                found = true;
                data = store;
                return;
            end

            part = parts{1};

            if ~isfield(store, part)
                found = false;
                data = [];
                return;
            end

            [found, data] = rooftopMoM.Cache.getNested_(store.(part), parts(2:end));
        end

        function store = clearNested_(store, parts)
            if isempty(parts)
                store = struct();
                return;
            end

            part = parts{1};

            if ~isfield(store, part)
                return;
            end

            if isscalar(parts)
                store = rmfield(store, part);
                return;
            end

            store.(part) = rooftopMoM.Cache.clearNested_(store.(part), parts(2:end));

            % If the child exists and became an empty struct, remove it from parent
            if isfield(store, part) && isstruct(store.(part)) && isempty(fieldnames(store.(part)))
                store = rmfield(store, part);
            end
        end
    end

    % Public static wrappers (call these directly)
    methods (Static)
        function setStatus(flag)
            inst = rooftopMoM.Cache.getInstance_();
            inst.setStatus_(flag);
        end

        function flag = getStatus()
            inst = rooftopMoM.Cache.getInstance_();
            flag = inst.getStatus_();
        end

        function setOption(path, flag)
            inst = rooftopMoM.Cache.getInstance_();
            inst.setOption_(path, flag);
        end

        function setOptionFromKey(key, flag)
            inst = rooftopMoM.Cache.getInstance_();
            inst.setOptionFromKey_(key, flag);
        end

        function flag = getOptionStatus(path)
            inst = rooftopMoM.Cache.getInstance_();
            flag = inst.getOptionStatus_(path);
        end

        function keys = listOptions()
            inst = rooftopMoM.Cache.getInstance_();
            keys = inst.listOptions_();
        end

        function autoSet(path, data)
            inst = rooftopMoM.Cache.getInstance_();
            inst.autoSet_(path, data);
        end

        function set(path, data)
            inst = rooftopMoM.Cache.getInstance_();
            inst.set_(path, data);
        end

        function data = get(path)
            inst = rooftopMoM.Cache.getInstance_();
            data = inst.get_(path);
        end

        function cleared = clear(path)
            inst = rooftopMoM.Cache.getInstance_();
            cleared = inst.clear_(path);
        end

        function snapshot = snapshot()
            inst = rooftopMoM.Cache.getInstance_();
            snapshot = inst.snapshot_();
        end

        function destroy()
            inst = rooftopMoM.Cache.getInstance_();
            inst.destroy_();
        end
    end
end
