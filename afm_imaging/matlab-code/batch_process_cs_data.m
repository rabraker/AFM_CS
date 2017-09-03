clear
clc
data_root = fullfile(getdataroot, 'cs-data');

fls = ls(fullfile(data_root, '*.csv'));

% ls results in one long string, separted by tab. Split that up. 
cs_file_list = strsplit(fls);

job_list = {};
purge_list = {};
for cs_data_file_cell = cs_file_list
    cs_data_file = cs_data_file_cell{1};
    
    cs_exp_meta_name = strrep(cs_data_file, '.csv', '-meta.mat');
    % Check if the meta file exists. If not, add it to a list of data to
    % delete.
    if exist(cs_exp_meta_name, 'file') ~= 2
        k = size(purge_list, 1);
        purge_list{k+1, 1} = cs_data_file;
        continue
    end

    % Check if the file has already been processed. If not, add it to the
    % job list. 
    img_data_file_name = strrep(cs_data_file, '.csv', '_img-data.mat');
    if exist(img_data_file_name, 'file') ~= 2
        k = size(job_list, 1);
        job_list{k+1,1} = cs_data_file;
    end
    
end

purge_list;

job_list; % 28

verbose = 0;

for job_file_cell = job_list'
    close all
    cs_job_file = job_file_cell{1};
    cs_exp_meta_name = strrep(cs_job_file, '.csv', '-meta.mat');
    
    % This is inefficient because we also load the file inside csdata2mat.
    % But I don't want to refactor right now.
    % Catch malformed .mat meta-data files.
    try
        load(cs_exp_meta_name)
        % Process the file:
        img_data  = csdata2mat(cs_job_file, cs_exp_meta_name, verbose);
    catch MExcp
        if strcmp(MExcp.identifier, 'MATLAB:textscan:EmptyFormatString')
            purge_list{length(purge_list)+1,1} = cs_job_file;
            fprintf('The file: \n %s \n is malformed. Adding parent to purge list.\n', cs_exp_meta_name)
            continue
        else
            rethrow(MExcp)
        end
        
    end
    
    % Construct the filename for image-data.
    img_data_file_name = strrep(cs_job_file, '.csv', '_img-data.mat');

    %Finally, save it.
    save(img_data_file_name, 'img_data')

   
    
end