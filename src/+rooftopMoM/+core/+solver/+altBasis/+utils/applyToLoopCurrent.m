function L = applyToLoopCurrent(L, loopCount, handle)
% applyToLoopCurrent: Returns the x vector whose part associated with loop currents is modified
%   by the user supplied handle.
%
%   Input:
%       L: Nx1 loop-tree current vector.
%       loopCount: Total number of loops on the antenna.
%       handle: Function handle that takes the loop part of the current vector.
%
%   Output:
%       L: Modified input loop current.

    idx = 1:loopCount;
    L(idx) = handle(L(idx));
end
