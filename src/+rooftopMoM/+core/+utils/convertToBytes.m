function bytes = convertToBytes(value, unit)
% Convert a numeric value + unit string into bytes.
%
%   Inputs:
%       value: Value to convert.
%       unit: Unit of the input value.
%
%   Output:
%       bytes: The inputed value expressed in bytes.

    unit = lower(unit);

    switch unit
        case {'b','bit'}
            bytes = value/8;
        case {'byte','bytes'}
            bytes = value;
        case {'kb','kilobyte'}
            bytes = value * 1e3;
        case {'kib','kibibyte'}
            bytes = value * 2^10;
        case {'mb','megabyte'}
            bytes = value * 1e6;
        case {'mib','mebibyte'}
            bytes = value * 2^20;
        case {'gb','gigabyte'}
            bytes = value * 1e9;
        case {'gib','gibibyte'}
            bytes = value * 2^30;
        case {'tb','terabyte'}
            bytes = value * 1e12;
        case {'tib','tebibyte'}
            bytes = value * 2^40;
        otherwise
            error('Unknown unit "%s".', unit);
    end
end
