clear;
clc;
% close all;


I = imread([pwd '/lena_256.jpg']);

if length(size(I)) > 2
I = rgb2gray(I);
end
I = double(I);


[n m] = size(I);
I_noisy = I + randn(n, n)*10;

E = muPathMaskGen(15,n,m,0.15);
% I = make_CS20NG(5, 5, 256, 10)*256
figure(1)
subplot(1,2,1)
imagesc(I_noisy)
colormap('gray')
%%
% clc
mu = 0.1*1;
lamy = 0.1/1000;
lamx = 0.1/50;
gamma = mu/10;
weight = 1; %0.001;%mu*5;
sigma = 1e-4;

N = 30;
tic
% [Ir] = BregmanSplitwithVerticalPenalty(I,E,0.03,0.0001,0.03,1000,40)
% [Ir] = BregmanSplitwithVerticalPenalty_delxy(I,E,weight,mu,lamy, lamx, gamma, 10, sigma);
pix_idx = find(CsTools.pixmat2vec(E) > 0.5);

[Ir] = bpvv_bregman(I_noisy,pix_idx,weight,mu,lamy, gamma, 10, sigma);
ttot = toc;
%%
figure(1)
subplot(1,2,2)
imagesc(Ir)
colormap('gray')

pval = max(max(I_noisy(:)), max(Ir(:)));
fprintf('total time: %.3f,  PSNR: %.3f,  SSIM:  %.3f\n', ttot,...
  psnr(Ir, I_noisy, pval), ssim(I_noisy, Ir, 'DynamicRange', pval))
%%
opts = l1qc_dct_opts('l1_tol', 0, 'epsilon', .1, 'cgmaxiter', 200);
pix_idx = find(E > 0.5);
b = I_noisy(pix_idx);

im_bp = l1qc_dct(256*256, 1, b, pix_idx, opts);
%
im_bp = reshape(im_bp, 256, []);
%%
mu = .1
im_bp_tv = SplitBregmanROF(im_bp, mu, 0.001);


figure(4)
subplot(1,2,1)
imagesc(im_bp)
colormap('gray')


subplot(1,2,2)
imagesc(im_bp_tv)
colormap('gray')

pval = max(max(I_noisy(:)), max(im_bp_tv(:)));
fprintf(' PSNR: %.3f,  SSIM:  %.3f\n', ...
  psnr(im_bp_tv, I_noisy, pval), ssim(I_noisy, im_bp_tv, 'DynamicRange', pval))
%%

%%
%     E_vec = PixelMatrixToVector(E);
%     pix_idx = find(E_vec==1);
%     I = speye(n^2);
%     R = I(pix_idx, :);
    
N = 256
sparsity = .25;
mu = .1/100;
lambda = mu;
gamma = mu/1000;

  % build an image of a square
% image = zeros(N,N);
% image(N/4:3*N/4,N/4:3*N/4)=255;
%  

 % build the sampling matrix, R
R = rand(N,N);
R = double(R<sparsity);

 % Form the CS data

F = R.*dct2(I);

% Recover the image
recovered = mrics(R,F, mu, lambda, gamma,30, 4);


imagesc(abs(recovered))
colormap('gray')


%%
I = eye(8);

clc
Dy(I)
Dyt(I)

function d = Dx(u)
[rows,cols] = size(u); 
d = zeros(rows,cols);
d(:,2:cols) = u(:,2:cols)-u(:,1:cols-1);
d(:,1) = u(:,1)-u(:,cols);
end

function d = Dxt(u)
[rows,cols] = size(u); 
d = zeros(rows,cols);
d(:,1:cols-1) = u(:,1:cols-1)-u(:,2:cols);
d(:,cols) = u(:,cols)-u(:,1);
end

function d = Dy(u)
[rows,cols] = size(u); 
d = zeros(rows,cols);
d(2:rows,:) = u(2:rows,:)-u(1:rows-1,:);
d(1,:) = u(1,:)-u(rows,:);
end







