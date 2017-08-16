clc
clear
pix = 256;
width = 5
pix_per_volt = (pix/width)*5

% dat = csvread('C:\Users\arnold\Documents\labview\afm_imaging/data/cs-data-out02.csv');
dat = csvread('C:\Users\arnold\Documents\labview\afm_imaging/data/data-out_ontable_csimage.csv');
ind_meas = find(dat(:, 6) ~= 0);

dat_meas = dat(ind_meas, :);
dat_meas(:,5) = detrend(dat_meas(:,5)); % detrend u_z as a long vector
% split each measurement into cells
%%
% bin all the data into pixels. 
for k = 1:max(dat_meas(:,end))
   inds = find(dat_meas(:,end) == k); 
   
   X_ks = dat_meas(inds, 1);
   xd = max(X_ks) - min(X_ks);
   npix = ceil(xd*pix_per_volt);
   
   Y_ks = dat_meas(inds, 2);
   E_ks = dat_meas(inds, 3);
   U_ks = dat_meas(inds, 5);
   
   x_pix = floor(pixelize(X_ks, npix)*256);
   y_pix = floor(pixelize(Y_ks, npix)*256);
   e_pix = pixelize(E_ks, npix);
   u_pix = pixelize(U_ks, npix);
  pix_dat{k} = [x_pix, y_pix, e_pix, u_pix];
end
%%
close all
figure; hold on;
hole_depth = (20/7)*(1/1000)*(20)
I = zeros(pix, pix);
pixelifsampled = I;
for k = 1:length(pix_dat)
    u_pix = pix_dat{k}(:,4);
    x_pix = pix_dat{k}(:,1);
    y_pix = pix_dat{k}(:,2);

    if abs(max(u_pix) - min(u_pix))> hole_depth*.5
        have_edge = 1;
        cs = 'r';
       % Then we have an edge. 
%         u_pix = u_pix - mean(u_pix);
        u_pix = u_pix - max(u_pix);
    else
        have_edge = 0;
        cs = 'b';
%         u_pix = u_pix - mean(u_pix);
        u_pix = u_pix - max(u_pix);
    end
    
%     figure;
    plot(u_pix, 'color', cs)
    ylim([-hole_depth, hole_depth])
    
%     keyboard
    
    for jj = 1:length(u_pix)
            n_row = y_pix(jj)+1;
            m_col = x_pix(jj)+1;
            I(n_row, m_col) = u_pix(jj);
            pixelifsampled(n_row, m_col) = 1;
    end
    
    
end
figure
lo = min(min(I));
hi = max(max(I));
imshow(flipud(I), [lo, hi])


%%

addpath('C:\Users\arnold\Documents\labview\reconstruction\BP\')
addpath('C:\Users\arnold\Documents\labview\reconstruction\BP\l1magic');
[n m] = size(I);

I_vector = PixelMatrixToVector(I);
pixelifsampled = PixelMatrixToVector(pixelifsampled);
I_vector = I_vector(find(pixelifsampled>0.5));

A = @(x) IDCTfun(x,pixelifsampled);
At = @(x) DCTfun(x,pixelifsampled);


close all;


Ir = idct(l1qc_logbarrier(At(I_vector), A, At, I_vector, 0.1));


close all;


Ir = real(Ir);
subplot(1,2,1)
imshow(I.*PixelVectorToMatrix(pixelifsampled,[n m]),[lo, hi]);
title('sample');
subplot(1,2,2)
imshow(PixelVectorToMatrix(Ir,[n m]),[min(Ir), max(Ir)]);
title('BP reconstruction');







