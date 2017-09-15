% function [sys_TF,num,den, sys_TF_delay] = fitDTTF2freqdat(FRF_data,freqs_rad,Ts,n,m,d,W)
% Least squares fit of freq data to Discrete-Time TF
%
% Author: Jeffrey Butterworth
% Date:   November 20, 2008
% n, m, d, % Order of: poles, zeros, delay
% FRF_data is complex data
%
% Modified: Arnold Braker
% Date:   March 25, 2016

function [sys_TF,num,den, sys_TF_delay] = fitDTTF2freqdat(FRF_data,freqs_rad,Ts,n,m,d,W)

%Weight constriction
W=diag([W;W]);


PHI = [-kron(ones(1,n),FRF_data).*exp(kron(-j*(1:n),Ts*freqs_rad)), ...
    exp(kron(-j*(d:d+m),Ts*freqs_rad));
           -kron(ones(1,n),conj(FRF_data)).*exp(kron(--j*(1:n),Ts*freqs_rad)), ...
    exp(kron(--j*(d:d+m),Ts*freqs_rad))];

% Recall... for PHI...
% We are considering the full spectrum of freqs, so we need:
% a) A portion of PHI that considers the positive freqs and standard FRF
% data AND
% b) A portion of PHI that considers the negative (two negs makes a pos)
% freqs and conj(FRF data)
% This helps us get REAl coeffs!
%===========================

% Least squares solution (short version)
% THETA_hat=(PHI)\([FRF_data; conj(FRF_data)]);


% THETA_hat=(W*PHI)\(W*[FRF_data; conj(FRF_data)]);
opts = optimset('lsqlin');
THETA_hat = lsqlin(W*PHI, W*[FRF_data; conj(FRF_data)], [], [], [], [], [], [], [], opts);

% THETA_hat=(PHI)\([FRF_data; conj(FRF_data)]);

den  = real([1; THETA_hat(1:n)].'); % Real is only needed to clear the near-zero...
num  = real([THETA_hat(n+1:end)].'); % ... imaginary terms.
den2 = real([den,zeros(1,m-n+d)].');

sys_TF = tf(num,den2.',Ts);
Nd     = length(zeros(1, m-n+d));
sys_TF_delay = tf(num, den, Ts, 'iodelay', Nd);

end




