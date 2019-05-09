
% close all
clear, clc
addpath('functions')
try; rmpath('functions/scanning_v0'); end
addpath('functions/scanning_v1')

rng(10)


npix = 512;  % image resolution.
sub_sample_frac = 0.15;  % Percent of pixels to subsample.

N_prescan = 250;
mu_pix_1 = 1
mu_pix_52 = 52




% ************************************************

[pix_mask_1, pix_idx_1] = mu_path_mask(mu_pix_1, npix, npix, sub_sample_frac, false);
[pix_mask_52, pix_idx_52] = mu_path_mask(mu_pix_52, npix, npix, sub_sample_frac, false);

%%
cs20ng = make_CS20NG(10, 10, npix, 10);

im = double(imread('standard_test_images/cameraman.tif'))/255;
% im = cs20ng;
I_1 = ones(npix,npix)-pix_mask_1;
I_52 = ones(npix,npix)-pix_mask_52;

%
F1 = mkfig(1, 9, 6); clf



ax = tight_subplot(2, 3, .01, .01, .01, true);
ax = reshape(ax', 3, 2)'

imagesc(ax(1,1), im), 
colormap('gray')

imagesc(ax(1,2), I_1)

colormap('gray')
remove_ticks(ax)

imagesc(ax(1,3), I_52)

colormap('gray')



cs_sim_1 = CsSim(im, pix_mask_1); 
cs_sim_52 = CsSim(im, pix_mask_52); 

cs_sim_1.solve_bp(true, true)
cs_sim_52.solve_bp(true, true)

%%
mu = 20
cla(ax(2,2))
img_1 = breg_anistropic_TV(cs_sim_1.Img_bp, mu, 0.001, 1000);
imagesc(ax(2,2), img_1);
colormap('gray')

img_52 = breg_anistropic_TV(cs_sim_52.Img_bp, mu, 0.001, 1000);
imagesc(ax(2,3), img_52);
colormap('gray')

remove_ticks(ax);

% save_fig(


ssim(img_1, im)
ssim(img_52, im)
%%
N = size(im, 1);
pix_idx = find(CsTools.pixmat2vec(pix_mask_52) > 0.5);
A = @(x) CsTools.Afun_dct(x, pix_idx);
At = @(x) CsTools.Atfun_dct(x, pix_idx, N, N);

% zr = zeros(N^2, 1);

M_fun = @(x) idct(x);
Mt_fun = @(x) dct(x);

E_fun = @(x) CsTools.E_fun1(x, pix_idx);
Et_fun = @(x) CsTools.Et_fun1(x, pix_idx, N, N);

b = CsTools.pixmat2vec(im);
b = b(pix_idx);


delta = 1e-2;
mu = 1e-2;

opts = NESTA_opts('Verbose', 10, 'errFcn', @(x)norm(x),...
    'U', M_fun, 'Ut', Mt_fun)
tic
[x,niter] = NESTA_mine(E_fun, Et_fun, b, mu, delta, opts);
toc

Im_est = CsTools.pixvec2mat(x, N);

ssim(Im_est, im)

imagesc(ax(2,1), Im_est);

%%
cs20_sorted = dct2(cs20ng);
cs20_sorted = sort(abs(cs20_sorted(:)), 'descend');

cam_sorted = dct2(im);
cam_sorted = sort(abs(cam_sorted(:)), 'descend');

figure(2); clf
semilogx(log10(cs20_sorted))
hold on
semilogx(log10(cam_sorted))


%%

skips = 1:4:512
if(skips(end) < npix)
    skips = [skips, npix];
end
img_sub = im(skips, :);

img_interp = zeros(npix, npix);
for k=1:npix

img_interp(:, k) = interp1(skips, img_sub(:, k), 1:512);


end


imagesc(ax(2,1), img_interp)
colormap('gray')



remove_ticks(ax)


ssim(img_interp, im)

