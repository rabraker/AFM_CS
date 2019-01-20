clear;
clc;
close all;


I = imread([pwd '/lena_256.jpg']);

if length(size(I)) > 2
I = rgb2gray(I);
end
I = double(I);


[n m] = size(I);


E = CsTools.muPathMaskGen(15,n,m,0.15);

%%
% clc
mu = 0.1*1;
lamy = 0.1/1000;
lamx = 0.1/500;
gamma = mu/10;
weight = 1; %0.001;%mu*5;
sigma = 1e-5;

% % mu = 0.1*10;
% % lambda = 0.1*1000;
% % gamma = mu/100;
% % weight = 0.00001;%mu*5;
% % sigma = 1e-5
N = 30;
tic
% [Ir] = BregmanSplitwithVerticalPenalty(I,E,0.03,0.0001,0.03,1000,40)
[Ir] = BregmanSplitwithVerticalPenalty_delxy(I,E,weight,mu,lamy, lamx, gamma, 10, sigma);
ttot = toc;
%%
pval = max(max(I(:)), max(Ir(:)));
fprintf('total time: %.3f,  PSNR: %.3f,  SSIM:  %.3f\n', ttot,...
  psnr(Ir, I, pval), ssim(I, Ir, 'DynamicRange', pval))

%%
figure(2)
imagesc(Ir)
colormap('gray')
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







