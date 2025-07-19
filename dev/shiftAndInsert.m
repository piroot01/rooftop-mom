% will be used in construction of the reordered impedance matrix in the rank 1 DOF add
    function idxList = shiftAndInsert(a, key)
% shiftAndInsert: Find spot for key in sorted a, then insert N+1 there
%   idxList = shiftAndInsert(a,key) returns a 1×(N+1) vector
%     [1,2,…,k, N+1, k+1,…,N]
%   in O(log N + N) time but with just two colon operations.

    N = numel(a);

    % Binary search to find k such that a(k) <= key < a(k+1) (or k=0 or k=N)
    lo = 1;  
    hi = N;
    while lo <= hi
        mid = floor((lo + hi)/2);
        if a(mid) <= key
            lo = mid + 1;
        else
            hi = mid - 1;
        end
    end
    k = hi;

    % Build 1:(N+1) once, then rotate the suffix in one go ——
    idxList = uint32(1):uint32(N+1);
    % shift everything from position k+1…end right by one, 
    % and drop the old last element into slot k+1
    idxList(k+1:end) = [N+1, idxList(k+1:end-1)];
end
