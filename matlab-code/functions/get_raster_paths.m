% [raster_paths] = get_raster_paths(data_root, data_name)
%
% Given data_root and data_name, returns a struct raster_paths with the 
% following fields:
% 
% raster_paths.data_path
% raster_paths.meta_path
% raster_paths.meta_in_path

function [raster_paths] = get_raster_paths(data_root, data_name)
  parent_root = get_parent_root(data_root);
  parent_name = get_parent_name(data_name, '_out_');
  
  raster_exp_meta_name = strrep(data_name, '.csv', '-meta.mat');
  raster_paths.data_path = fullfile(data_root, data_name);
  raster_paths.meta_path = fullfile(data_root, raster_exp_meta_name);
  
  raster_paths.parent_path = fullfile(parent_root, parent_name);
  k = regexp(raster_paths.data_path, '_out');
  raster_paths.meta_in_path = sprintf('%s.mat', raster_paths.data_path(1:k-1));
end