% will be used in construction of the reordered impedance matrix in the rank 1 DOF add

function idxList = shiftAndInsert(a, key)
    N = numel(a);
    lo = 1;
    hi = N;
    % Binary search using bit-shift for faster midpoint calculation
    while lo <= hi
        mid = bitshift(lo + hi, -1); % Equivalent to floor((lo+hi)/2) but faster
        if a(mid) <= key
            lo = mid + 1;
        else
            hi = mid - 1;
        end
    end
    k = hi; % k is in [0, N]
    
    % Build the index list in a vectorized way (conditionals removed)
    idxList = zeros(1, N+1);
    idxList(1:k) = 1:k;        % Safe even if k=0 (no-op)
    idxList(k+1) = N+1;        % Insert the new element
    idxList(k+2:end) = k+1:N;  % Safe even if k=N (no-op)
end

