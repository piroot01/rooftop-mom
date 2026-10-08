function [output, applicationTime] = measureApplicationTime(applicationHandle, repetitionCount)
% measureApplicationTime: Measures the mean time of a repeated operator application.
%
%   Inputs:
%       applicationHandle: Function handle without inputs returning the application output.
%       repetitionCount: Positive integer specifying the number of the timed repetitions.
%
%   Outputs:
%       output: Output of the application.
%       applicationTime: Mean time of a single application in seconds.

    % Warm-up application
    output = applicationHandle();

    tic;
    for n = 1:repetitionCount
        output = applicationHandle();
    end
    applicationTime = toc / repetitionCount;
end
