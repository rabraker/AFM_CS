clear
Tfig_path = @(fname) PATHS.tuesday_fig_path('12-11-2018', fname);


% Task 1: DCT of full sim image. Threshold and reconstruct with usual pix mask. 
% Look up the results on K-sparsity. Does it make sense to make it k-sparse in
% the sense that if k= 0.1 * 512^2?

% What is the difference between thresholding the 1D-dct and the 2D-dct??
% clear
npix = 256
x_start = 26/2; 
y_start = 26/2;
img_mat = make_CS20NG(x_start, y_start, npix);

figure(1)
imshow(img_mat, [0,1])
%%
clc

mu_pix = 1;


E = CsTools.muPathMaskGen(15,npix, npix,0.15);
pix_mask = CsTools.pixmat2vec(E);
% pix_mask = mu_path_mask(npix, mu_pix, 0.12);
pix_idx = find(pix_mask == 1);

figure(2)
imagesc(E)
colormap('gray')
sum(pix_mask(:))/npix^2

At = @(x) CsTools.Atfun_dct(x,pix_idx, npix^2); %E^T*M^T

xorig = CsTools.pixmat2vec(img_mat);
%%
clc
a = .7*1.8;
t = (0:1:npix^2-1)';

g = tf(a, [1 a]);
xx = lsim(g, xorig, t);
xx = lsim(g, flipud(xx), t);
xx = flipud(xx);
figure(10); clf
plot(t, xorig)
hold on
plot(t, xx, '--')

ylim([-.05, 1.05])
% xlim([2.54*0, 2.545]*1e5)
xlim([3250, 3550])
%%
b = xx(pix_idx);
x0 = At(b);

tic;
opts = CsTools.l1qc_opts('warm_start_cg', 0);
eta_rec = CsTools.l1qc(x0, b, pix_idx-1, opts);
t_l1qc = toc();
fprintf('mex time: %.3f\n', t_l1qc)
%%
% Xorig = img_mat;
Xorig = CsTools.pixvec2mat(xx, npix);

x_rec = idct(eta_rec)';

f1 = @()SB_ITV(x_rec(:), .5);
f2 = @() SplitBregmanROF(CsTools.pixvec2mat(x_rec, npix),5,0.001);

timeit(f1)
timeit(f2)

%%


x_rec_filt = SB_ITV(x_rec(:), .5);
X_l1qc_filt = SplitBregmanROF(CsTools.pixvec2mat(x_rec, npix),5,0.001);

X_l1qc_filt = reshape(X_l1qc_filt,npix,npix);
X_l1qc = CsTools.pixvec2mat(x_rec, npix);
% X_l1qc_filt = CsTools.pixvec2mat(x_rec_filt, npix);


dnr = max(max(img_mat(:)), max(X_l1qc(:)));
dnr_f = max(max(img_mat(:)), max(X_l1qc_filt(:)));
psnr_l1qc = psnr(Xorig, X_l1qc, dnr)
psnr_l1qc_f = psnr(Xorig, X_l1qc_filt, dnr_f)
ssim_l1qc = ssim(Xorig, X_l1qc, 'DynamicRange', dnr)
ssim_l1qc_f = ssim(Xorig, X_l1qc_filt, 'DynamicRange', dnr_f)

fprintf('time l1qc: %5f,  l1qc psnr: %.3f l1qc ssim: %.3f \n', t_l1qc, psnr_l1qc, ssim_l1qc)
fprintf('time l1qc: %5f,  l1qc psnr: %.3f l1qc ssim: %.3f \n', t_l1qc, psnr_l1qc_f, ssim_l1qc_f)
figure(3)
subplot(1,3,1)
imagesc(X_l1qc)
colormap('gray')

subplot(1,3,2)
imagesc(X_l1qc_filt)
colormap('gray')


subplot(1,3,3)
imagesc(Xorig)
colormap('gray')

%%

% % For BregmanSplit
% mu = 0.1;
% gamma = 1*01;   %mu/100;
% lam = 1*mu;
% epsy = 0.0011; %0.001;%mu*5;
% epsx = 0.0011;
% sigma = 1e-5;

% For BregmanSplit_delxy
mu = 0.2;
gamma = .1*01;   %mu/100;
lam = 1*mu;
epsy = 01; %0.001;%mu*5;
epsx = 01;
sigma = 1e-5;

tic
% [Ir] = BregmanSplitwithVerticalPenalty(I,E,0.03,0.0001,0.03,1000,40)
% [Ir] = BregmanSplit(Xorig, E, mu, gamma, 10, sigma);
[Ir] = CsTools.BregmanSplit_delxy(Xorig, E, epsy, epsx, mu, lam, gamma, 10, sigma);
ttot = toc;

pval = max(max(img_mat(:)), max(Ir(:)));

fprintf('\n\n------------------------------------------------------\n')
fprintf('time l1qc: %5f,  l1qc psnr: %.3f l1qc ssim: %.3f \n', t_l1qc, psnr_l1qc, ssim_l1qc)
fprintf('time (SB): %.3f,  SB PSNR: %.3f   SB SSIM:  %.3f\n', ttot,...
  psnr(Ir, img_mat, pval), ssim(img_mat, Ir, 'DynamicRange', pval))