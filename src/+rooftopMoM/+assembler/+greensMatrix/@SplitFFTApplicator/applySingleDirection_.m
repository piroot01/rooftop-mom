function output = applySingleDirection_(~, input, evenBlock, oddBlock, modulation, demodulation)
% applySingleDirection_: Applies the Green's matrix of a strip, the rooftops of one direction only.
%
%   Inputs:
%       input: Nx1 vector of the rooftop coefficients.
%       evenBlock: Spectral block of the even subband (a column for an x strip, a row for a y strip).
%       oddBlock: Spectral block of the odd subband, shaped as evenBlock.
%       modulation: Odd subband modulation, shaped as evenBlock.
%       demodulation: Complex conjugate of the modulation.
%
%   Output:
%       output: Nx1 vector of the product.

    count = numel(input);

    % Zero-extend the input to the subband size and orient it as the blocks.
    input(numel(evenBlock)) = 0;
    input = reshape(input, size(evenBlock));

    output = ifft(evenBlock .* fft(input));

    % The odd subband is transformed in the input buffer, which is not needed afterwards.
    input = ifft(oddBlock .* fft(input .* modulation)) .* demodulation;
    output = output + input;

    output = reshape(output(1:count), [], 1);
end
