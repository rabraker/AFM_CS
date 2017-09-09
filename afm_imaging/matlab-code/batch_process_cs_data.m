clear
clc
addpath('functions')
data_root = fullfile(getdataroot, 'cs-data');
fls = ls(fullfile(data_root, '*.csv'));

% The ls command works differently on linux vs windows. What the fuck
% MATLAB??? How hard could it have been to make the output the same accross
% platforms??
% Account for that.
if ispc
    cs_file_list = {};
    for k = 1:size(fls,1) % m by n char array
        cs_file_list{k} = fullfile(data_root, fls(k,:));
    end
else
    % ls results in one long string, separted by tab. Split that up. 
    cs_file_list = strsplit(fls);
end



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
%%

for job_file_cell = job_list'
    close all
    cs_job_file = job_file_cell{1};
    cs_exp_meta_name = strrep(cs_job_file, '.csv', '-meta.mat');
    
    k = regexp(cs_job_file, 'Hz_out')
    
    meta_in_path = sprintf('%s.mat', cs_job_file(1:k+1));
    % This is inefficient because we also load the file inside csdata2mat.
    % But I don't want to refactor right now.
    % Catch malformed .mat meta-data files.
    try
        load(cs_exp_meta_name)
        if exist(meta_in_path, 'file') == 2
            load(meta_in_path)
        else % all new data has a meta.mat file. all old data had same params we need
            CsExpMetaIn = struct('npix', 256, 'width', 5)
            keyboard
        end
        % Process the file:
        img_data  = csdata2mat(cs_job_file, cs_exp_meta_name,CsExpMetaIn, verbose);
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



for fl_cell = cs_file_list
    cs_exp_data_path = fl_cell{1};
    
    mat_data_path = strrep(cs_exp_data_path, '.csv', '_img-data.mat');
    
    
    if exist(mat_data_path, 'file') ~= 2
        continue
    end
    if exist(meta_in_path, 'file') ~= 2
        continue
    end
    
%     fig_root = 'C:\Users\arnold\Documents\labview\afm_imaging\matlab-code\figures';
    fig_root = getfigfigroot();
    cs_exp_fig_name = strrep(cs_exp_data_path, '.csv', '-fig.fig');
   
    if exist(cs_exp_fig_name, 'file') == 2
        continue
    end
    
    load(mat_data_path)
    
    width = img_data.width;
    Ts = img_data.Ts;
    pix = 256;
    
    f5 = figure(6); clf
    subplot(2,3,1)
    ax3 = gca();
    imshow_sane(img_data.cs_im, ax3, width, width);
    title('sample');

    bp_im = detrend_plane(img_data.bp_im);
    subplot(2,3,2)
    ax4 = gca();
    % imshow_sane(PixelVectorToMatrix(Ir_bp,[n m]), ax4, width, width);
    imshow_sane(img_data.bp_im, ax4, width, width);
    title('BP reconstruction');

    subplot(2,3,3)
    ax5 = gca();
    imshow_sane(img_data.smp_im, ax5, width, width)
    title('SMP reconstruction');

    fname = strsplit(cs_exp_fig_name, filesep);
    fname = fname{end};
    s = metadata2text(img_data.meta, Ts);
    s1 = sprintf('%s\nperc=%.3f', s, sum(sum(img_data.pixelifsampled))/pix/pix);
    s2 = sprintf('file: %s\n', fname);
%         disp(s1);
%         disp(s2);
    
        subplot(2,3, [4,5,6]);
        ax4 = gca();
        ax4.Visible = 'off';
        t1 = text(0,.5, s1, 'Units', 'normalized');
        t2 = text(0, t1.Extent(2)-.1, s2, 'Units', 'normalized', 'interpreter', 'none');


    saveas(f5, cs_exp_fig_name, 'fig')
    
end




