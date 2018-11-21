function img_data = get_acc_cs_img(perc)
  
  data_root = '/media/labserver/acc2018-data/';
% dat_root = '/media/labserver/acc2018-data/data/';


% -------- Constants --------
Ts = 40e-6;
cs_root =  fullfile(data_root,'cs-data', '5microns', '9-22-2017');
if perc == 5
  cs_name = 'cs-traj-512pix-5perc-500nm-5mic-01Hz_out_9-23-2017-04.csv';
elseif perc == 10
  cs_name = 'cs-traj-512pix-10perc-500nm-5mic-01Hz_out_9-23-2017-04.csv';
elseif perc == 15
  cs_name = 'cs-traj-512pix-15perc-500nm-5mic-01Hz_out_9-23-2017-03.csv';
else
  error('Must supply perc = (5|10|15) (as double)');
end


% width = 20;
% xdata = [0, width];
% ydata = [0, width];

% Prep file names
img_data_name = strrep(cs_name, '.csv', '_img-data.mat');
parent_name = get_parent_name(cs_name, '_out_');
meta_name = strrep(parent_name, '.csv', '.mat');

% load meta-data and data. 
load(fullfile(cs_root, meta_name)); % supplies CsExpMetaIn

load(fullfile(cs_root, img_data_name)); % supplies img_data
img_data.CsExpMetaIn = CsExpMetaIn;
pixelifsampled = img_data.CsExpMetaIn.pixelifsampled;
pix = img_data.CsExpMetaIn.npix;
actual_sub_sample_frac = length(find(pixelifsampled == 1))/pix^2;


img_data.perc = actual_sub_sample_frac*100;
img_data.bp_im_normalized = img_data.bp_im - mean(img_data.bp_im(:));

img_data.thresh = (20/7)*(1/1000)*22;
  
  
end