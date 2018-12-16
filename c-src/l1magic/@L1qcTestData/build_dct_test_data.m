function build_dct_test_data(test_data_root)
  % Create a small set of test data

  EMx_path = fullfile(test_data_root, 'dct_small_EMx.json');
  MtEty_path = fullfile(test_data_root, 'dct_small_MtEty.json');
  MtEt_Emx_path = fullfile(test_data_root, 'dct_small_MtEt_EMx.json');
  N = 50;
  Ts = 1/50;
  To = 0.25;
  
  omega =2*pi/To;
  
  k = [ 0:N-1]';
  
  x = sin(omega * Ts *k);
  
  pix_idx = [2, 10, 15, 20, 25, 30, 35, 40, 45, 49];
  
  Ny = length(pix_idx);
  pix_mask = zeros(N,1);
  pix_mask(pix_idx) = 1;
  
  yy = IDCTfun(x, pix_mask);
  
  data = struct('x0', x(:)', 'x1', yy(:)', 'pix_idx', pix_idx(:)'-1);
  savejson('', data, EMx_path);
  
  
  xx = DCTfun(yy, pix_mask);
  xx_save = xx;
   
  data = struct('x0', xx_save(:)', 'x1', yy(:)', 'pix_idx', pix_idx(:)'-1);
  savejson('', data, MtEty_path);
  %
  data = struct('x0', x(:)', 'x1', xx_save(:)', 'pix_idx', pix_idx(:)'-1);
  savejson('', data, MtEt_Emx_path);
  
  %%
  
  % Now, the real deal:
  fpath = '/home/arnold/matlab/afm-cs/matlab-code/notes/data/cs_sim_CS20NG.mat';
  addpath ~/matlab/afm-cs/reconstruction/BP

  dat = load(fpath);
  cs_sim = dat.cs_sim;

  pix_mask_vec = PixelMatrixToVector(cs_sim.pix_mask);
  
  
  
  pix_idx = find(pix_mask_vec > 0.5);
  
  % y = E*M*x
  img_vec = PixelMatrixToVector(cs_sim.Img_sub_sampled);
  % y, set of measurements. have to remove all the spots we didn't sample.
  y_vec = img_vec(pix_idx);
  A = @(x) IDCTfun(x, pix_mask_vec); % E*M
  At = @(y) DCTfun(y, pix_mask_vec); %E^T*M^T

  EMx = A(img_vec);
  MtEty = At(y_vec);
  
  MtEt_EMx = At(A(img_vec));
  
  jopts.FloatFormat = '%.15f';
  jopts.FileName = fullfile(test_data_root, 'dct_large.json');
  savejson('', struct('x_in', img_vec(:)', 'y_in', y_vec(:)',...
    'EMx', EMx(:)', 'MtEty', MtEty(:)', 'MtEt_EMx', MtEt_EMx(:)', 'pix_idx', pix_idx(:)'-1), jopts);
  
  
  
  
  
% For debugging, just the dct part. 
%   fpath = '/home/arnold/matlab/afm-cs/matlab-code/notes/data/cs_sim_CS20NG.mat';
%   addpath ~/matlab/afm-cs/reconstruction/BP
% 
%   dat = load(fpath);
%   cs_sim = dat.cs_sim;
% 
%   img_vec = PixelMatrixToVector(cs_sim.Img_sub_sampled);
% %   img_vec = img_vec(1:2048);
%   pix_idx = [1:length(img_vec)];
%   
%   % y, set of measurements. have to remove all the spots we didn't sample.
%   y_vec = img_vec;
%   A = @(x) idct(x);
%   At = @(y) dct(y);
% 
%   EMx = A(img_vec);
%   MtEty = At(y_vec);
%   
%   MtEt_EMx = At(A(img_vec));
%   
%   jopts.FloatFormat = '%.20f';
%   jopts.FileName = fullfile(test_data_root, 'dct_large.json');
%   savejson('', struct('x_in', img_vec(:)', 'y_in', y_vec(:)',...
%     'EMx', EMx(:)', 'MtEty', MtEty(:)', 'MtEt_EMx', MtEt_EMx(:)', 'pix_idx', pix_idx(:)'-1), jopts);
%     
  

end