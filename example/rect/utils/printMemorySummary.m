function printMemorySummary(names, methods, storedMemory, peakMemory)
% printMemorySummary: Prints a compact ASCII table comparing the memory of the operator applications.
%
%   Inputs:
%       names: 1xM cell of the strategy names.
%       methods: 1xM cell of the method descriptions.
%       storedMemory: 1xM vector of the memory in bytes of the stored operators.
%       peakMemory: 1xM vector of the peak working memory in bytes of one application (output included).

    line = '+--------------------+----------------------------------+-----------------+-------------------+';
    fprintf('%s\n', line);
    fprintf('| %-18s | %-32s | %-15s | %-17s |\n', 'Strategy', 'Method', 'Stored (MB)', 'Peak working (MB)');
    fprintf('%s\n', line);

    for n = 1:numel(names)
        fprintf('| %-18s | %-32s | %15.3f | %17.3f |\n', ...
            names{n}, methods{n}, storedMemory(n) / 2^20, peakMemory(n) / 2^20);
    end

    fprintf('%s\n', line);
end
