

function [cs_paths] = get_cs_paths(data_root, data_name)
  
  parent_root = get_parent_root(data_root);
  parent_name = get_parent_name(data_name, '_out_');
  parent_meta_name = strrep(parent_name, '.csv', '.json');
  
  cs_exp_meta_name = strrep(data_name, '.csv', '-meta.json');
  cs_paths.data_path = fullfile(data_root, data_name);
  cs_paths.parent_path = []; %fullfile(parent_root, parent_name);
  cs_paths.parent_meta_path = fullfile(parent_root, parent_meta_name);
  cs_paths.meta_path = fullfile(data_root, cs_exp_meta_name);
 % k = regexp(cs_paths.data_path, '_out');
  cs_paths.meta_in_path = []; %sprintf('%s.mat', cs_paths.data_path(1:k-1));
  
  cs_paths.data_path_mat = [] ;% strrep(cs_paths.data_path, '.csv', '.json');
end