classdef CacheKeys
% CacheKeys: Central registry of cachable variable names (human-readable constants).
%
%   Usage examples:
%       % dot-path constant:
%       keyDot = rooftopMoM.CacheKeys.CORE_IMPEDANCE_FULLASSEMBLER_ROOFTOPPOINTS;
%
%       % convert to underscore-style base (without _enable/_disable suffix)
%       base = rooftopMoM.CacheKeys.toUnderscore(keyDot); % 'core_impedanceMatrix_FullAssembler_rooftopPoints'
%
%       % full enable/disable keys:
%       enableKey = rooftopMoM.CacheKeys.enableKeyFromDot(keyDot); % 'core_impedanceMatrix_FullAssembler_rooftopPoints_enable'
%       disableKey = rooftopMoM.CacheKeys.disableKeyFromDot(keyDot); 'core_impedanceMatrix_FullAssembler_rooftopPoints_disable'
%
%       % list all defined dot-paths:
%       allKeys = rooftopMoM.CacheKeys.listAll();

    properties (Constant)
        % Naming convention: UPPER_SNAKE with readable descriptive name.
        CORE_IMPEDANCEMATRIX_FULLASSEMBLER_QUADRATUREPOINTS = 'core.impedanceMatrix.FullAssembler.quadraturePoints';
        CORE_IMPEDANCEMATRIX_FULLASSEMBLER_ROOFTOPPOINTS = 'core.impedanceMatrix.FullAssembler.rooftopPoints';
        CORE_IMPEDANCEMATRIX_FULLASSEMBLER_ROOFTOPVALUES = 'core.impedanceMatrix.FullAssembler.rooftopValues';

        CORE_GREENSMATRIX_FULLASSEMBLER_QUADRATUREPOINTS = 'core.greensMatrix.FullAssembler.quadraturePoints';
        CORE_GREENSMATRIX_FULLASSEMBLER_ROOFTOPPOINTS = 'core.greensMatrix.FullAssembler.rooftopPoints';
        CORE_GREENSMATRIX_FULLASSEMBLER_ROOFTOPVALUES = 'core.greensMatrix.FullAssembler.rooftopValues';
    end

    methods (Static)
        function underscorePath = toUnderscore(dotPath)
        % Convert a dot-path -> underscore base (no _enable/_disable)

            underscorePath = strrep(char(dotPath), '.', '_');
        end

        function key = enableKeyFromDot(dotPath)
        % Return underscore-style enable key for setOptionFromKey

            base = rooftopMoM.CacheKeys.toUnderscore(dotPath);
            key = [base, '_enable'];
        end

        function key = disableKeyFromDot(dotPath)
        % Return underscore-style disable key for setOptionFromKey

            base = rooftopMoM.CacheKeys.toUnderscore(dotPath);
            key = [base, '_disable'];
        end

        function dotPath = dotFromUnderscore(underscoreBase)
        % Convert underscore base back to dot-path (useful if user provides underscore)

            dotPath = strrep(char(underscoreBase), '_', '.');
        end

        function s = listAll()
        % Return a struct mapping constant names -> dot-path values.

            mc = metaclass(rooftopMoM.CacheKeys);
            props = mc.PropertyList;
            s = struct();

            for k = 1:numel(props)
                p = props(k);
                if p.Constant && strcmp(p.GetAccess, 'public')
                    name = p.Name;
                    s.(name) = rooftopMoM.CacheKeys.(name);
                end
            end
        end
    end
end
