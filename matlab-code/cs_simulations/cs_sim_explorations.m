% This script looks at different CS schemes for the CS20-NG simulation grating.
% In particular, I look at using the standard 52pixel mu-path length for
% several sampling percentages. Then I compare that to the ideal of 1-pixel
% length. 
%
% Also, I look what happens when we threshold the master image in the DCT
% domain.
%
% And finally, what happens when we do the same CS schemes based on the
% 10%thresholded DCT as "master"
%

Tfig_path = @(fname) PATHS.tuesday_fig_path('12-11-2018', fname);


%% Task 1: DCT of full sim image. Threshold and reconstruct with usual pix mask. 
% Look up the results on K-sparsity. Does it make sense to make it k-sparse in
% the sense that if k= 0.1 * 512^2?

% What is the difference between thresholding the 1D-dct and the 2D-dct??
% clear
x_start = 26; 
y_start = 26;
img_mat = make_CS20NG(x_start, y_start);

imshow(img_mat, [0,1])

%% Task 2: With same threshold, use 4 pix masks
% -1. Fully random (mulength = 1)
% -2. Mu length = .25 of usual
%   ''  .5
%   ''  .75  of usual

mu_pix = 52;
frac_s = [0.05, 0.12, 0.25, .35];
cs_sims = cell(1, length(frac_s));
parfor iter = 1:length(frac_s)
  pix_mask = mu_path_mask(npix, mu_pix, frac_s(iter));
  
  cs_sim_iter = CsSim(img_mat, pix_mask);
  cs_sim_iter.solve_bp;
  cs_sims{iter} = cs_sim_iter;
end

% NOW SET, MU-PIX TO 1 PIXEL.
iter = 1;
mu_pix = 1;
%%
cs_sims_2 = cell(1, length(new_fracs));
parfor iter = 1:length(frac_s)
%   pix_mask = mu_path_mask(npix, mu_pix, frac_s(iter));

  pix_mask = rand(npix, npix);

  idx_keep = find(pix_mask < new_fracs(iter));
  idx_toss = find(pix_mask >= new_fracs(iter));
  pix_mask(idx_keep) = 1;
  pix_mask(idx_toss) = 0;

  cs_sim_iter = CsSim(img_mat, pix_mask);
  cs_sim_iter.solve_bp;
  cs_sims_2{iter} = cs_sim_iter;
end


%%
% -- Now, plot everything.
f22 = figure(22); clf
[ha] = tight_subplot(2, 4, [.03, .005 ], .027, .005);

new_fracs = frac_s*0; %
for k=1:length(frac_s)
  mu_pix = 52;
  ax_iter = ha(k);

  imagesc(ax_iter, cs_sims{k}.Img_bp, [0, 1]);
  colormap('gray')
  set(ax_iter, 'YTick', [], 'XTick', [])
  
  sample_frac = sum(cs_sims{k}.pix_mask(:))/npix^2;
  new_fracs(k) = sample_frac;
  ssim_iter = ssim(img_mat, cs_sims{k}.Img_bp);
  psnr_iter = psnr(img_mat, cs_sims{k}.Img_bp);
  
  stit = sprintf('mu-pix = %d, $%.1f\\%%$, SSIM=%.3f, PSNR=%.1f', mu_pix,...
    sample_frac*100, ssim_iter, psnr_iter);
  title(ax_iter, stit)
end



for k=1:length(frac_s)
  mu_pix = 1;
  ax_iter = ha(4 + k);
    
  imagesc(ax_iter, cs_sims_2{k}.Img_bp, [0, 1]);
  colormap('gray')
  set(ax_iter, 'YTick', [], 'XTick', [])
  
  ssim_iter = ssim(img_mat, cs_sims_2{k}.Img_bp);
  psnr_iter = psnr(img_mat, cs_sims_2{k}.Img_bp);
  
  sample_frac = sum(cs_sims_2{k}.pix_mask(:))/npix^2;
  stit = sprintf('mu-pix = %d, $%.1f\\%%$, SSIM=%.3f, PSNR=%.1f', mu_pix,...
    sample_frac*100, ssim_iter, psnr_iter);
  
  title(ax_iter, stit)
end
%%
saveas(f22, Tfig_path('cs_sim_mu_path_vs_singlepix.fig'));
save_fig(f22, PATHS.note_fig('cs_sim_mu_path_vs_singlepix'));

%%
img_vec = PixelMatrixToVector(img_mat);
img_vec_dct = dct(img_vec);

npix = 512;
perc_s = [1, 0.19, 0.096, 0.05, 0.01, 0.001]
K_s = perc_s*0;
img_mats_thresh = cell(length(perc_s),1);
for iter = 1:length(perc_s)
K = floor(npix^2 * perc_s(iter));
K_s(iter) = K;
img_vec_thresh = idct(k_sparsify(img_vec_dct, K));

img_mat_thresh = PixelVectorToMatrix(img_vec_thresh, [npix, npix]);
img_mats_thresh{iter} = img_mat_thresh;
% figure
% imshow(img_mat_thresh, [0, 1])
end

%
f2 = figure(2); clf
[ha, pos] = tight_subplot(2, 3, [.03, .005 ], .027, .005);


f3 = figure(3);
figure(3); clf
ax3 = gca();
plot(ax3, log10(img_vec_dct_mag_sorted));
hold(ax3, 'on');

h_s = gobjects(length(perc_s),1);
for k=1:length(perc_s)
  ax_iter = ha(k);
  imagesc(ax_iter, img_mats_thresh{k}, [0, 1])
  colormap('gray')
  set(ax_iter, 'YTick', [], 'XTick', [])
  thresh_frac = perc_s(k);
  stit = sprintf('threshhold perc=%.1f', thresh_frac*100);
  title(ax_iter, stit)
  
  h_s(k) = plot([K_s(k),K_s(k)], ylim, '--');
  h_s(k).DisplayName = stit;
  
  
end
legend(h_s)
title('sorted DCT coeficcients (log10)')

%%
saveas(f2, Tfig_path('cs20ng_sim_thresholds_img.fig'));
save_fig(f2, PATHS.note_fig('cs20ng_sim_thresholds_img'));


saveas(f3, Tfig_path(cs20ng_sim_thresholds_mag.fig'));
save_fig(f3, PATHS.note_fig('cs20ng_sim_thresholds_mag'));


%%
% ---------------------------------------------------------------------------- %
% ---------- No, do the same experiment, but with the thresholded image. ----- %
% % Im going to take the 10% guy.
ten_perc_idx = 3;
img_mat_thresh10 = img_mats_thresh{ten_perc_idx};
mu_pix = 52;
frac_s = [0.05, 0.12, 0.25, .35];
cs_sims_thresh_52pix = cell(1, length(frac_s));
parfor iter = 1:length(frac_s)
  pix_mask = mu_path_mask(npix, mu_pix, frac_s(iter));
  
  cs_sim_iter = CsSim(img_mat_thresh10, pix_mask);
  cs_sim_iter.solve_bp;
  cs_sims_thresh_52pix{iter} = cs_sim_iter;
end

% NOW SET, MU-PIX TO 1 PIXEL.
iter = 1;
mu_pix = 1;

cs_sims_thresh_1pix = cell(1, length(new_fracs));
parfor iter = 1:length(frac_s)
%   pix_mask = mu_path_mask(npix, mu_pix, frac_s(iter));

  pix_mask = rand(npix, npix);

  idx_keep = find(pix_mask < new_fracs(iter));
  idx_toss = find(pix_mask >= new_fracs(iter));
  pix_mask(idx_keep) = 1;
  pix_mask(idx_toss) = 0;

  cs_sim_iter = CsSim(img_mat_thresh10, pix_mask);
  cs_sim_iter.solve_bp;
  cs_sims_thresh_1pix{iter} = cs_sim_iter;
end


%%
% -- Now, plot everything.
f32 = figure(32); clf
[ha, pos] = tight_subplot(2, 4, [.03, .005 ], .027, .005);

new_fracs = frac_s*0; %
for k=1:length(frac_s)
  mu_pix = 52;
  ax_iter = ha(k);

  imagesc(ax_iter, cs_sims_thresh_52pix{k}.Img_bp, [0, 1]);
  colormap('gray')
  set(ax_iter, 'YTick', [], 'XTick', [])
  
  sample_frac = sum(cs_sims_thresh_52pix{k}.pix_mask(:))/npix^2;
  new_fracs(k) = sample_frac;

  ssim_iter = ssim(img_mat_thresh10, cs_sims_thresh_52pix{k}.Img_bp);
  psnr_iter = psnr(img_mat_thresh10, cs_sims_thresh_52pix{k}.Img_bp);
   stit = sprintf('mu-pix = %d, $%.1f\\%%$, SSIM=%.3f, PSNR=%.1f', mu_pix,...
    sample_frac*100, ssim_iter, psnr_iter);
  
  title(ax_iter, stit)
end



for k=1:length(frac_s)
  mu_pix = 1;
  ax_iter = ha(4 + k);
    
  imagesc(ax_iter, cs_sims_thresh_1pix{k}.Img_bp, [0, 1]);
  colormap('gray')
  set(ax_iter, 'YTick', [], 'XTick', [])
  
  
  sample_frac = sum(cs_sims_thresh_1pix{k}.pix_mask(:))/npix^2;
  
  ssim_iter = ssim(img_mat_thresh10, cs_sims_thresh_1pix{k}.Img_bp);
  psnr_iter = psnr(img_mat_thresh10, cs_sims_thresh_1pix{k}.Img_bp);
  stit = sprintf('mu-pix = %d, $%.1f\\%%$, SSIM=%.3f, PSNR=%.1f', mu_pix,...
    sample_frac*100, ssim_iter, psnr_iter);
  
  title(ax_iter, stit)
end
%%
saveas(f32, Tfig_path('cs_sim_mu_path_vs_singlepix_thresh.fig'))
save_fig(f32, PATHS.note_fig('cs_sim_mu_path_vs_singlepix_thresh'));




function x_sparse = k_sparsify(x, k)
  [~,idx_descend] = sort(abs(x), 'descend');
  z = abs(x(idx_descend));
  z_thresh = z(k);
  
  x(abs(x)<z_thresh) = 0;
  x_sparse = x;
end


% % Now, lets try the same thing with the boat image.
% img_mat_boat = double(imread('boat.png')); % This is also 512 x 512
% npix = size(img_mat_boat,1);
% 
% K = floor(npix^2 * 0.025);
% img_vec_boat = PixelMatrixToVector(img_mat_boat);
% 
% img_boat_vec_dct = dct(img_vec_boat);
% img_boat_vec_dct_mag_sorted = sort(abs(img_boat_vec_dct), 'descend');
% 
% img_vec_boat_thresh = idct(k_sparsify(img_boat_vec_dct, K));
% 
% img_boat_mat_thresh = PixelVectorToMatrix(img_vec_boat_thresh, [npix, npix]);
% 
% 
% figure(2)
% ax1 = subplot(3,1,[1,2]);
% ax2 = subplot(3,1,3);
% imshow_dataview(img_boat_mat_thresh, [0, 255], ax1, ax2)
% 
% figure(3)
% hold on
% plot(log10(img_boat_vec_dct_mag_sorted/255));

% [LoD,HiD,LoR,HiR] = wfilters('bior3.5');
% [c1, s1] = wavedec(img_vec_boat, 5, LoD, HiD);
% 
% figure(3)
% plot(log10(sort(abs(c1)/255, 'descend')));
% %
% ylim([-10, 5])
% ax = gca();
% ax.XScale = 'log';