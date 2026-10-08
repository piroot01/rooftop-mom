clc; clearvars;

N = 100000000;

% input preparation
ht = complex(randn([N, 1]), randn([N, 1]));
x = complex(randn([N, 1]), randn([N, 1]));

%% 1D reference multiplication
%t = toeplitz(ht, ht);

%y_ref = t * x;

%% 1D fft accelerated muliplication
function y = FFTMult_(ht, x)
    N = numel(ht);

    hc = [ht; 0; flipud(ht(2:end))];
    
    y = ifft(fft(hc) .* fft(x, 2*N));
    y = y(1:N);
end

y_fftMult = FFTMult_(ht, x);

%% 1D split fft accelerated multiplication
function y = splitFFTMult_(ht, x)
    N = numel(ht);
    
    hc = [ht; 0; flipud(ht(2:end))];

    hc = fft(hc);
    hce = hc(1:2:end);
    hco = hc(2:2:end);

    k = (0:(N - 1)).';
    P = exp(-1i * pi * k ./ N);

    ye = ifft(hce .* fft(x));
    yo = ifft(hco .* fft(P .* x)) .* conj(P);

    y = 0.5 * (ye + yo);
end

y_splitFFTMult = splitFFTMult_(ht, x);

%% 1D compare the reference and fft accuracy
%diff_ref_FFTMult = norm(y_ref - y_fftMult)
%diff_ref_splitFFTMult = norm(y_ref - y_splitFFTMult)
diff_FFTMult_splitFFtMult = norm(y_fftMult - y_splitFFTMult)