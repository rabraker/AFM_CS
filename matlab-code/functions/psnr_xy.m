% [psnr] = psnr_xy(X, Y)
%
% Implement peak signal to noise ratio, from 
% A comparison of reconstruction methods for undersampled atomic force
% microscopy images, Yufan Luo and Sean B Andersson, Nanotech. 2015
  

function [psnr] = psnr_xy(X, Y)

[n1, n2] = size(X);

sum_diff = 0;
for i_row = 1:n1
   for j_col = 1:n2
       sum_diff = sum_diff + (X(i_row, j_col) - Y(i_row, j_col))^2;
   end
end

psnr = sqrt((1/n1/n2)*sum_diff);
psnr = max(max(X))/psnr;
psnr = 20*log10(psnr);




end

