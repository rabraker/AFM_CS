
data_root = fullfile(getdataroot, 'cs-data');
data_in = 'cs-traj-10perc-1000nm-10mic-01Hz.csv';
exist(fullfile(data_root, data_in), 'file')
%%
% create meta file name
meta_in = strrep(data_in, '.csv', '.mat');
meta_data_path = fullfile(data_root, meta_in);

width = 10;
nom_perc = 10;
raster_freq = 1;
mu_length = 1000;
tip_velocity = width/(0.5*(1/raster_freq));

CsExpMetaIn.width = width;
CsExpMetaIn.nom_perc = nom_perc;
CsExpMetaIn.mu_length = mu_length;
CsExpMetaIn.tip_velocity = tip_velocity;
CsExpMetaIn.npix = 256;

save(meta_data_path, 'CsExpMetaIn')