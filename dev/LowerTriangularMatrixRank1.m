n = 10000;
L = tril(rand(n));
x = rand(n,1);

A = rand(n, n);

tic;
y1 = ltmv(L, x);
t1 = toc;
tic;
y2 = L*x;
t2 = toc;
fprintf('ltmv time: %.4f s, full multiply time: %.4f s\n', t1, t2);

function y = ltmv(L, x)
    % Validate inputs
    [~, n] = size(L);

    % Preallocate output
    y = zeros(n,1);

    % Main loop: forward substitution style
    for i = 1:n
        % Only multiply up to the diagonal
        y(i) = L(i,1:i) * x(1:i);
    end
end
