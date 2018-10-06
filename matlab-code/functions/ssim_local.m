% [ ssim] = structural_sim(x, y)
%
% Implements the Structural Similarity index from 
% Image quality assessment: from error visibility to structural similarity,
% Zhou Wang and A. C. Bovik and H. R. Sheikh and E. P. Simoncelli, IEEE 
% Trans. Image Processing, 2004

function [ ssim_xy] = ssim_local(X, Y)

x = PixelMatrixToVector(X);
y = PixelMatrixToVector(Y);

% from paper. 
K1 = 0.01;
K2 = 0.03;
L = 255;

mu_x = mean(x);
mu_y = mean(y);

sig_x = sqrt(var(x));
sig_y = sqrt(var(y));

cov_xy = cov(x,y);
sig_xy = cov_xy(1,2);
% keyboard

L = 255;

C1 = (K1*L)^2;
C2 = (K2*L)^2;
C3 = C2/2;

l_xy = (2*mu_x*mu_y + C1)/(mu_x^2 + mu_y^2 + C1);
c_xy = (2*sig_x*sig_y + C2)/(sig_x^2 + sig_y^2 + C2);
s_xy = (sig_xy + C3)/(sig_x*sig_y + C3);

% from paper
alpha = 1;
beta = 1;
gamma = 1;

ssim_xy = (l_xy^alpha)*(c_xy^beta)*(s_xy^gamma);

end

