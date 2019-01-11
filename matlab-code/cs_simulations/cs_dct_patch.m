clear
Tfig_path = @(fname) PATHS.tuesday_fig_path('12-11-2018', fname);


% Task 1: DCT of full sim image. Threshold and reconstruct with usual pix mask. 
% Look up the results on K-sparsity. Does it make sense to make it k-sparse in
% the sense that if k= 0.1 * 512^2?

% What is the difference between thresholding the 1D-dct and the 2D-dct??
% clear
x_start = 26; 
y_start = 26;
img_mat = make_CS20NG(x_start, y_start);

imshow(img_mat, [0,1])
%%
clc
mu_pix = 1;
npix = size(img_mat, 1);

patch_sz = 64;

M_patches = dct_M_patches(npix, patch_sz);

img_dct_patch = dct2_patched(img_mat, M_patches, 2);

figure(2)
imagesc(img_dct_patch);
colormap('gray')

sp_dct = k_sparsify(img_dct_patch, floor(0.1*512^2));
imag_again = dct2_patched(sp_dct, M_patches, 3);

figure(3)
imagesc(imag_again)
colormap('gray')

pix_mask = mu_path_mask(npix, mu_pix, 0.12);
pix_idx = find(pix_mask == 1);

figure(4)
imagesc(pix_mask)
colormap('gray')
sum(pix_mask(:))/npix^2


A_ = @(x) A_fun(x, M_patches, pix_idx, npix);
At_ = @(x) At_fun(x, M_patches, pix_idx, npix);
b = PixelMatrixToVector(img_mat);
b = b(pix_idx);
x0 = At_(b);

eta_rec = CsTools.l1qc_logbarrier(x0, A_, At_, b, 0.1);
eta_mat = PixelVectorToMatrix(real(eta_rec), [npix, npix]);
%%
X_rec = dct2_patched(eta_mat, M_patches, 3);

figure(4)
imagesc(X_rec)
colormap('gray')

function x = At_fun (y, M_patches, pix_idx, N)
  
  Ety = zeros(N*N,1);
  Ety(pix_idx) = y;
  
  ETY = CsTools.pixvec2mat(Ety, N);
  
  X = dct2_patched(ETY, M_patches, 2);
  x = CsTools.pixmat2vec(X);
  
end

function y = A_fun (x, M_patches, pix_idx, N)
  X = CsTools.pixvec2mat(x, N);
  Y = dct2_patched(X, M_patches, 3);
  y = CsTools.pixmat2vec(Y);
  y = y(pix_idx);
  
end

function M = dct_M_patches(N, patch_sz)
  N_patch = N/patch_sz;
  M_blk = sparse(dct_M(patch_sz));
  I = sparse(eye(N_patch));
  M = kron(I, M_blk);
end

function M = dct_M(N)
  
  E = eye(N);
  M = dct(E, [], 1);

end

function img_dct_patch = dct2_patched(img_mat, M_patches, xfm)

  if xfm ==2
    img_dct_patch = M_patches * img_mat * M_patches';
  elseif xfm ==3
    img_dct_patch = M_patches' * img_mat * M_patches;
  end

end
% function img_dct_patch = dct2_patched(img_mat, patch_sz, xfm)
%   npix = size(img_mat,1);
%   img_dct_patch = img_mat*0;
%   n_patch = npix/patch_sz;
%   assert (n_patch - floor(n_patch) == 0);
% 
%   if xfm ==2
%     dct_fun = @(x) dct2(x);
%   elseif xfm ==3
%     dct_fun = @(x) idct2(x);
%   end
%   for col=1:n_patch
%     idx_c = (col-1)*patch_sz + 1;
%     for row = 1:n_patch
%       idx_r = (row-1)*patch_sz + 1;
% 
%       patch_dct = dct_fun(img_mat(idx_r:idx_r+patch_sz-1, idx_c:idx_c+patch_sz-1));
% 
%       img_dct_patch(idx_r:idx_r+patch_sz-1, idx_c:idx_c+patch_sz-1) = patch_dct;
%     end
%   end
%   
% end
function x_sparse = k_sparsify(x, k)
  [~,idx_descend] = sort(abs(x), 'descend');
  z = abs(x(idx_descend));
  z_thresh = z(k);
  
  x(abs(x)<z_thresh) = 0;
  x_sparse = x;
end
