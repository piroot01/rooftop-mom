function bytes = generateByteCount(varargin)
% generateByteCount: Returns the bytes occupied by the given arrays.
%
%   Input:
%       varargin: Arrays of any class.
%
%   Output:
%       bytes: Total number of bytes.

    bytes = 0;

    for n = 1:numel(varargin)
        array = varargin{n}; %#ok<NASGU> Read by whos.
        information = whos('array');
        bytes = bytes + information.bytes;
    end
end
