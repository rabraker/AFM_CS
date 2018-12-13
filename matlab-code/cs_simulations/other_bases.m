clear
fpath = '/home/arnold/matlab/afm-cs/matlab-code/notes/data/cs_sim_CS20NG.mat';
addpath ~/matlab/afm-cs/reconstruction/BP
addpath(genpath('~/matlab/dependencies/SparseLab2.1-Core/Utilities/'))
load(fpath)
%%
clc
clc
npix = 512;
frac = .01;
img_mat = cs_sim.Img_original;
pix_mask_vec = PixelMatrixToVector(cs_sim.pix_mask);
Idx = find(pix_mask_vec > 0.5);

qmf = MakeONFilter('Haar', 5);
L = 2
inv_fun = @(x) IWT_PO(x, L, qmf);
fwd_fun = @(x) FWT_PO(x, L, qmf);

inv_fun2 = @(x) IWT2_PO(x, L, qmf);
fwd_fun2 = @(x) FWT2_PO(x, L, qmf);

% ------ Compare images with truncated sparsity in different bases. 
img_vec = PixelMatrixToVector(cs_sim.Img_original);
pix_sparse_dct_vec = sparsify_1d_signal(img_vec, frac*npix^2, @dct, @idct);
img_dct1_sparse = PixelVectorToMatrix(pix_sparse_dct_vec, [npix, npix]);

pix_sparse_fht1_vec = sparsify_1d_signal(img_vec, frac*npix^2, fwd_fun, inv_fun);
img_fht1_sparse = PixelVectorToMatrix(pix_sparse_fht1_vec, [npix, npix]);

img_dct2_sparse = sparsify_2d_signal(img_mat, frac*npix^2, @dct2, @idct2);

img_fht2_sparse = sparsify_2d_signal(img_mat, frac*npix^2, fwd_fun2, inv_fun2);


figure(1); clf
% subplot(2,2, 1)
% imagesc(cs_sim.Img_original, [0, 0.06])
% colormap('gray')
subplot(2,2, 1)
imagesc(img_dct1_sparse, [0, 0.06])
ssim_1 = ssim(img_dct1_sparse, img_mat);
psnr_1 = psnr(img_dct1_sparse, img_mat);
colormap('gray')
title(sprintf('DCT-1D, PSNR=%.3f, SSIM=%.3f', psnr_1, ssim_1))

subplot(2,2, 2)
imagesc(img_fht1_sparse, [0, 0.06])
ssim_2 = ssim(img_fht1_sparse, img_mat);
psnr_2 = psnr(img_fht1_sparse, img_mat);
colormap('gray')
title(sprintf('FHT-1D, PSNR=%.3f, SSIM=%.3f', psnr_2, ssim_2))

subplot(2,2, 3)
imagesc(img_dct2_sparse, [0, 0.06])
ssim_3 = ssim(img_dct2_sparse, img_mat);
psnr_3 = psnr(img_dct2_sparse, img_mat);
colormap('gray')
title(sprintf('dct-2D, PSNR=%.3f, SSIM=%.3f', psnr_3, ssim_3))


subplot(2,2, 4)
imagesc(img_fht2_sparse, [0, 0.06])
ssim_4 = ssim(img_fht2_sparse, img_mat);
psnr_4 = psnr(img_fht2_sparse, img_mat);
colormap('gray')
title(sprintf('FHT-2D, PSNR=%.3f, SSIM=%.3f', psnr_4, ssim_4))

%%
% inv_fun = @(x)idct(x)
% fwd_fun = @(x)dct(x);
% A = @(x)Afun(x, Idx, inv_fun);
% At = @(y)Atfun(y, Idx, npix^2, fwd_fun);
% 
% img_bp_other = cs_sim.solve_bp_other(A, At);
% img_vec = idct(img_bp_other);
% img = PixelVectorToMatrix(img_vec, [npix, npix]);
% figure(2)
% imshow(img, [0,0.06])
%%
inv_fun = @(x)FHT2(x);
fwd_fun = @(x)FHT2(x);

A = @(x)Afun(x, Idx, inv_fun);
At = @(y)Atfun(y, Idx, npix^2, fwd_fun);

img_bp_other = cs_sim.solve_bp_other(A, At);

img_vec = inv_fun(img_bp_other);
img = PixelVectorToMatrix(img_vec, [npix, npix]);
figure(2)
imshow(img, [0,0.06])

img_vec_filt = TIDenoise(img_vec', 'S', qmf);
img_filt = PixelVectorToMatrix(img_vec_filt, [npix, npix]);

figure(3)
imshow(img_filt, [0,0.06])


%%
% qmf = MakeONFilter('Haar',1);
% L = floor(log2(log(npix^2)))+1;

qmf = MakeONFilter('Daubechies', 8);
% L = floor(log2(log(npix^2)))+1;
% L = 2^4
L = 2
inv_fun = @(x) IWT_PO(x, L, qmf);
fwd_fun = @(x) FWT_PO(x, L, qmf);


A = @(x)Afun(x, Idx, inv_fun);
At = @(y)Atfun(y, Idx, npix^2, fwd_fun);

img_bp_other = cs_sim.solve_bp_other(A, At);
%%
img_vec = inv_fun(img_bp_other);
img = PixelVectorToMatrix(img_vec, [npix, npix]);
figure(5)
imshow(img, [0,0.06])

img_vec_filt = TIDenoise(img_vec', 'S', qmf);
img_filt = PixelVectorToMatrix(img_vec_filt, [npix, npix]);

figure(6)
imshow(img_filt, [0,0.06])
%%
frac = 0.1

pix_mask_vec = rand(npix^2, 1);

idx_keep = find(pix_mask_vec < frac);
idx_toss = find(pix_mask_vec >= frac);
pix_mask_vec(idx_keep) = 1;
pix_mask_vec(idx_toss) = 0;
pix_mask = PixelVectorToMatrix(pix_mask_vec, [npix, npix]);

cs_sim_iter = CsSim(cs_sim.Img_original, pix_mask);
Idx = find(pix_mask_vec > 0.5);

% qmf = MakeONFilter('Symmlet', 8);
qmf = MakeONFilter('Daubechies', 8);
L = 2;
inv_fun = @(x) IWT2_PO(x, L, qmf);
fwd_fun = @(x) FWT2_PO(x, L, qmf);

% inv_fun = @(x) idct2(x);
% fwd_fun = @(x) dct2(x);

A = @(x)Afun2(x, Idx, npix, inv_fun);
At = @(y)Atfun2(y, Idx, npix, fwd_fun);

img_bp_other = cs_sim_iter.solve_bp_other(A, At);
%%
img_vec = inv_fun(PixelVectorToMatrix(img_bp_other, [npix, npix]));
img = PixelVectorToMatrix(img_vec, [npix, npix]);
figure(5)
ax1 = subplot(3,1,[1,2]);
ax2 = subplot(3,1,3);
imshow_dataview(img, [0,0.06], ax1, ax2)
%%
img_vec_filt = TIDenoise(PixelMatrixToVector(img_vec)', 'H');
img_filt = PixelVectorToMatrix(img_vec_filt, [npix, npix]);

figure(6)
ax1 = subplot(3,1,[1,2]);
ax2 = subplot(3,1,3);
imshow_dataview(img_filt, [0,0.06], ax1, ax2)

function y = Afun2(x, Idx, n, inv_fun)
  Xmat = PixelVectorToMatrix(x, [n,n]);
  y = PixelMatrixToVector(inv_fun(Xmat));
  y = y(Idx);
end
function x = Atfun2(y, Idx, n, fwd_fun)
  x = zeros(n^2,1);
  
  x(Idx) = y;
  Xmat = PixelVectorToMatrix(x, [n,n]);
  x = PixelMatrixToVector(fwd_fun(Xmat));
end


function y = Afun(x, Idx, inv_fun)
  y = inv_fun(x);
  y = y(Idx);
end
function x = Atfun(y, Idx, n, fwd_fun)
  x = zeros(n,1);
  x(Idx) = y;
  x = fwd_fun(x);
end




