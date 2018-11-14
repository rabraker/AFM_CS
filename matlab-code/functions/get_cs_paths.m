

function [cs_paths] = cs_exp_paths(data_root, data_name)
  
  cs_exp_meta_name = strrep(data_name, '.csv', '-meta.mat');
  cs_paths.data_path = fullfile(data_root, data_name);
  cs_paths.meta_path = fullfile(data_root, cs_exp_meta_name);
  
  k = regexp(cs_paths.data_path, '_out');
  cs_paths.meta_in_path = sprintf('%s.mat', cs_paths.data_path(1:k-1));
end