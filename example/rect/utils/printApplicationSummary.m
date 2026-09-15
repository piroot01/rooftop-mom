function printApplicationSummary(names, methods, times, errors)
% printApplicationSummary: Prints a compact ASCII table comparing the operator applications. The first row is the
%   reference for the speed-up.
%
%   Inputs:
%       names: 1xM cell of the strategy names.
%       methods: 1xM cell of the method descriptions.
%       times: 1xM vector of the mean application times in seconds.
%       errors: 1xM vector of the relative errors with respect to the reference application.

    line = '+--------------------+----------------------------------+------------+------------+------------+';
    fprintf('%s\n', line);
    fprintf('| %-18s | %-32s | %-10s | %-10s | %-10s |\n', 'Strategy', 'Method', 'Time (ms)', 'Speed-up', 'Rel. error');
    fprintf('%s\n', line);

    for n = 1:numel(names)
        fprintf('| %-18s | %-32s | %10.4f | %10.2f | %10.2e |\n', ...
            names{n}, methods{n}, 1e3 * times(n), times(1) / times(n), errors(n));
    end

    fprintf('%s\n', line);
end
