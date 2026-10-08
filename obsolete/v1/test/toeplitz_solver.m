% example of a fast toep solver

T = [1, 2, 3; 2, 1, 2; 3, 2, 1];
e = 1;

Tinv = inv(T);

Tinv(:, 1)

invToeplitzColBGY(T(:, 1), e)

function x = invToeplitzColBGY(firstCol, j)
% invToeplitzColBGY  Fast O(n log n) solver for the j-th column of T^{-1}
%   using the Brent–Gustavson–Yun divide-and-conquer algorithm specialized
%   for unit-vector right-hand sides (monomials).
%
%   x = invToeplitzColBGY(firstCol, j)
%   Input:
%     firstCol - n-by-1 vector, first column of Toeplitz T
%     j        - integer in [1..n], column index of inverse desired
%   Output:
%     x        - n-by-1 vector satisfying T*x = e_j

    n = length(firstCol);
    if j<1 || j>n || floor(j)~=j
        error('j must be integer between 1 and n');
    end
    % Build firstRow
    firstRow = [firstCol(1); flipud(firstCol(2:end))];
    % Call recursive BGY solver
    x = bgy_solve(firstCol, firstRow, j-1);
end

function x = bgy_solve(c, r, shift)
% Recursive core: solve Toeplitz(c,r)*x = z^shift
    n = length(c);
    if n <= 32
        % Base case: form small T and solve directly
        T = toeplitz(c, r);
        b = zeros(n,1); b(shift+1) = 1;
        x = T\b;
        return;
    end
    % Split size
    m = floor(n/2);
    % Extract top-left m-by-m block params
    cA = c(1:m);
    rA = r(1:m);
    % Extract bottom-right (n-m)-by-(n-m) block params
    cB = c(m+1:-1:m+1-(n-m)+1);
    rB = r(m+1:-1:m+1-(n-m)+1);
    cB = flipud(cB);
    rB = flipud(rB);
    % Obtain symbol polynomials A(z), B(z)
    A = cA(:)'; B = cB(:)';
    % Fast Extended Euclidean on A and B: U*A + V*B = 1
    [U, V] = poly_ext_gcd(A, B);
    % Determine shifts into subproblems
    % RHS for first half: -V(z)*z^shift mod A, but shift then mod by truncating
    Vshift = shift_coeffs(V, shift, m);
    % RHS index j0 in first half if nonzero
    j0 = find(Vshift,1)-1;
    if isempty(j0), j0 = -1; end
    x0 = zeros(m,1);
    if j0>=0
        x0 = bgy_solve(cA, rA, j0);
        x0 = -x0;
    end
    % RHS for second half: U(z)*z^shift mod B
    Ushift = shift_coeffs(U, shift, n-m);
    j1 = find(Ushift,1)-1;
    x1 = zeros(n-m,1);
    if j1>=0
        x1 = bgy_solve(cB, rB, j1);
    end
    % Merge x0 and x1 into full solution
    x = [x0; x1];
end

function s = shift_coeffs(p, k, L)
% Shift polynomial p(z) by z^k and truncate/pad to length L
    s = zeros(1,L);
    for i=0:length(p)-1
        idx = i + k + 1;
        if idx>=1 && idx<=L
            s(idx) = p(i+1);
        end
    end
    s = s(:)';
end

function [U, V] = poly_ext_gcd(A, B)
% Extended GCD of polynomials A and B over R
% Returns U,V such that U*A + V*B = 1 (assuming A,B coprime)
    if all(B==0)
        U = 1./A(1);
        V = zeros(size(B));
        return;
    end
    % Euclid loop with simple polynomial remainder
    U0 = 1; V0 = zeros(size(B));
    U1 = zeros(size(A)); V1 = 1;
    a = A; b = B;
    while any(b)
        [q, r] = deconv(a, b);
        r = r(find(r,1):end);
        U2 = U0 - conv(q, U1);
        V2 = V0 - conv(q, V1);
        a = b; b = r;
        U0 = U1; U1 = U2;
        V0 = V1; V1 = V2;
    end
    % Normalize leading coefficient to 1
    lc = a(find(a,1));
    U = U0./lc; V = V0./lc;
end

