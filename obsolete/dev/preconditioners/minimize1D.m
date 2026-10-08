function [xmin, fmin, history] = minimize1D(f, x0, options)
% MINIMIZE1D 1-D minimizer (golden-section or greedy).
%
% [xmin, fmin, history] = minimize1D(f, x0, options)
%
% Inputs:
%   f       - function handle: scalar = f(scalar)
%   x0      - initial guess (scalar)
%   options - struct (all fields optional):
%             .method     = 'golden' (default) or 'greedy'
%             .tol        = 1e-6 (stopping tolerance on x interval or step)
%             .maxIter    = 200
%             .lb, .ub    = scalar bounds (if given, golden uses them; if not given, auto-bracket)
%             .init_step  = 0.1 (used for bracketing/greedy)
%             .step_red   = 0.5  (greedy)
%             .verbose    = false
%
% Outputs:
%   xmin    - found x (approximate minimizer)
%   fmin    - f(xmin)
%   history - struct with fields:
%               .x : visited x-values (column)
%               .f : corresponding f-values
%
% Examples:
%   f = @(x) (x-2).^2 + sin(5*x).*0.1;
%   [xm,fm,h] = minimize1D(f, 0, struct('method','golden','verbose',true));
%   [xm,fm,h] = minimize1D(f, 0, struct('method','greedy','init_step',0.5,'verbose',true));
%

    if nargin < 3, options = struct(); end
    method    = getOpt(options, 'method', 'golden');
    tol       = getOpt(options, 'tol', 1e-6);
    maxIter   = getOpt(options, 'maxIter', 200);
    lb        = getOpt(options, 'lb', -inf);
    ub        = getOpt(options, 'ub',  inf);
    init_step = getOpt(options, 'init_step', 0.1);
    step_red  = getOpt(options, 'step_red', 0.5);
    verbose   = getOpt(options, 'verbose', false);

    x0 = x0(:); % ensure scalar-like
    if numel(x0) ~= 1
        error('x0 must be scalar for 1-D minimization.');
    end
    x0 = double(x0);

    % wrapper to guard against NaN/Inf: treat them as +Inf
    safef = @(x) safe_eval(f, x);

    history.x = [];
    history.f = [];

    switch lower(method)
        case 'golden'
            % if bounds provided, use them; otherwise attempt to bracket
            if isfinite(lb) && isfinite(ub)
                a = lb; b = ub;
                if a >= b, error('lb must be < ub'); end
            else
                % attempt to find a bracket [a,b] with a < b and a minimum inside
                [a, b] = bracket_search(safef, x0, init_step, maxIter);
                if isempty(a)
                    error('Unable to bracket a minimum from x0. Try giving lb/ub or different x0/init_step.');
                end
            end
            if verbose, fprintf('Bracket: [%.6g, %.6g]\n', a, b); end

            % golden-section
            [xmin, fmin, hx, hf] = golden_section_search(safef, a, b, tol, maxIter, verbose);
            history.x = hx(:);
            history.f = hf(:);

        case 'greedy'
            % simple 1-D greedy hill-climb with step reduction
            x = max(min(x0, ub), lb);
            step = init_step;
            fcur = safef(x);
            history.x = x;
            history.f = fcur;
            iter = 0;
            while iter < maxIter && step > tol
                iter = iter + 1;
                improved = false;

                for signv = [1, -1]
                    cand = x + signv * step;
                    % enforce bounds
                    cand = max(min(cand, ub), lb);
                    fcand = safef(cand);
                    if fcand < fcur
                        x = cand;
                        fcur = fcand;
                        improved = true;
                        history.x(end+1,1) = x;
                        history.f(end+1,1) = fcur;
                        if verbose, fprintf('Iter %d: moved to x=%.6g f=%.6g (step=%.3g)\n', iter, x, fcur, step); end
                        break; % try from new x
                    end
                end

                if ~improved
                    step = step * step_red;
                    if verbose, fprintf('Iter %d: no improvement, reduce step -> %.3g\n', iter, step); end
                end
            end
            xmin = x; fmin = fcur;

        otherwise
            error('Unknown method "%s". Use "golden" or "greedy".', method);
    end
end

%% ---------------- Helper functions ----------------

function val = getOpt(s, name, def)
    if isfield(s, name), val = s.(name); else val = def; end
end

function y = safe_eval(f, x)
    try
        y = f(double(x));
        if ~isfinite(y) || ~isscalar(y)
            y = Inf;
        end
    catch
        y = Inf;
    end
end

function [a, b] = bracket_search(f, x0, step, maxIter)
% BRACKET_SEARCH expand from x0 to find a,b such that f has a minimum in [a,b]
% We look for a triple where f(x0) is lower than neighbors or expand until we find
    fa = f(x0);
    % check neighbors
    x1 = x0 - step; x2 = x0 + step;
    f1 = f(x1); f2 = f(x2);
    if fa <= f1 && fa <= f2
        % already a local minimizer, return small bracket
        a = x1; b = x2;
        return;
    end

    % choose expansion direction: toward the side that decreases
    if f1 < fa
        dir = -1;
    else
        dir = 1;
    end

    a = x0; fa = f(x0);
    b = x0 + dir * step; fb = f(b);
    iter = 0;
    while iter < maxIter
        iter = iter + 1;
        % expand exponentially
        step = step * 2;
        c = b + dir * step;
        fc = f(c);
        % want bracket a < b < c with f(b) < f(a) and f(b) < f(c)
        if dir > 0
            a = a; b = b; % ensure a < b < c
            if fa >= fb && fb <= fc
                % success: ensure ordering a < b < c
                a = x0; b = c; % we'll use [a,c] as bracket that contains optimum near old b
                % more conservative: return [x0-step, c] maybe fine
                a = x0 - abs(step); b = c;
                return;
            end
        else
            % dir < 0
            if fc >= fb && fb <= fa
                % we found a minimum at b, return bracket [c, x0+abs(step)]
                a = c; b = x0 + abs(step);
                return;
            end
        end
        % shift window
        a = b;
        fa = fb;
        b = c;
        fb = fc;
    end
    % failed
    a = []; b = [];
end

function [xopt, fopt, hx, hf] = golden_section_search(f, a, b, tol, maxIter, verbose)
% GOLDEN_SECTION_SEARCH bracket [a,b] known to contain a unimodal minimum.
    phi = (sqrt(5)-1)/2;        % ~0.618...
    c = b - phi*(b-a);
    d = a + phi*(b-a);
    fc = f(c);
    fd = f(d);
    hx = [c; d];
    hf = [fc; fd];

    iter = 0;
    while (b-a) > tol && iter < maxIter
        iter = iter + 1;
        if fc < fd
            b = d;
            d = c;
            fd = fc;
            c = b - phi*(b-a);
            fc = f(c);
            hx(end+1,1) = c; hf(end+1,1) = fc;
        else
            a = c;
            c = d;
            fc = fd;
            d = a + phi*(b-a);
            fd = f(d);
            hx(end+1,1) = d; hf(end+1,1) = fd;
        end
        if verbose && mod(iter,10)==0
            fprintf('Gold it %d: interval = [%.6g, %.6g]\n', iter, a, b);
        end
    end

    % pick best among sampled points
    [fopt, idx] = min(hf);
    xopt = hx(idx);
end