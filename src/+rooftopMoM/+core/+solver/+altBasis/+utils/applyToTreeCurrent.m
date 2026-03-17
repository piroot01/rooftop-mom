function T = applyToTreeCurrent(T, loopCount, handle)
% applyToTreeCurrent: Returns the x vector whose part associated with tree currents is modified
%   by the user supplied handle.
%
%   Input:
%       T: Nx1 loop-tree current vector.
%       loopCount: Total number of loops on the antenna.
%       handle: Function handle that takes the tree part of the current vector.
%
%   Output:
%       T: Modified input tree current.

    idx = (loopCount + 1):length(T);
    T(idx) = handle(T(idx));
end
