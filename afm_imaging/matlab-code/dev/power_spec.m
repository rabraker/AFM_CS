

function [X, f_s] = power_spec(x, Ts)

    L = length(x);
    X_fft = fft(x);

    P2 = abs(X_fft/L);
    X = P2(1:L/2+1);
    X(2:end-1) = 2*X(2:end-1);

    f_s = (1/Ts)*[0:(L/2)]'/L;

    % df = pwelch(sin(50000*x_ref.time));
%     semilogx(f, 20*log10(X))



end