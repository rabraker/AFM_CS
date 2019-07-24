%% Generate x-y input waveforms.
%
% Use 0.2 hz x-dir triangle wave.
clear
clc
addpath('functions')
addpath(fullfile(getMatPath(), 'dependencies', 'jsonlab'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Where to save raster data

addpath('functions/state_space_x')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xdirControl = get_xdir_loop_shaped_control(false);
%%
% Trace is one line at sec_line. The whole period is trace and re-trace.
image_side = 5; % micro-meters.
target_dir = sprintf('%dmicrons', image_side);
data_root = PATHS.raster_image_data('5microns', 'parents-loop');

rast_mode = 1;
if rast_mode == 1
    rast_mode_name = '';
    % data_root = PATHS.raster_image_data('5microns', 'parents-loop');
        data_root = PATHS.raster_image_data('5microns', 'parents-loop-subL');
elseif rast_mode == 2
    rast_mode_name = '-subL';
    data_root = PATHS.raster_image_data('5microns', 'parents-loop-subL');
end



lines = [512, 128, 64];
freqs = [1, 2.5, 5, 8, 10];
% lines = 128;
% freqs = 1;
for pp = 1:length(lines)
%     npix = 128;
    npix = lines(pp);
    for k=1:length(freqs)
        % raster_freq =
        raster_freq = coerce_raster_freq(freqs(k)); % Hz.
        

        if rast_mode == 1
            rast = RasterTraj(image_side, raster_freq, npix);
        elseif rast_mode == 2
            rast = RasterTraj(image_side, raster_freq, npix, 'y_traj_gen', @subl_y_TG);
        end
        
        x_traj = rast.x_traj_volts;
        y_traj = rast.y_traj_volts;
        
        n_harmonics = 100;
        w_max = 170 * 2 * pi;
        [r_ff, r_des] = fourier_tri(x_traj.Data, raster_freq*2*pi, n_harmonics, xdirControl.Hy_rprime, w_max);
        
        [y, t] = lsim(xdirControl.Hyr, r_ff, x_traj.Time);
        
        figure(1); clf; hold on
        h1 = plot(x_traj.Time, x_traj.Data);
        h2 = plot(t, y, '--b');
        
        h4 = plot(t, r_ff, '-g');
        h3 = plot(t, r_des, '--m');
        
        h1.DisplayName = 'triangle';
        h2.DisplayName = 'stage-output';
        h3.DisplayName = 'smoothed triangle';
        h4.DisplayName = 'smoothed tri + Ginv';
        legend([h1, h2, h3, h4])
        
        % We have to interleave the x & y data like
        % [x(1), y(1), x(2), y(2), ....]
        % So, reshape
        % x = [1, 2, 3, 4, 5;
        %      6, 7, 8, 9, 10];
        xy_data = [x_traj.Data(:)';
                   y_traj.Data(:)'];
        xy_data = reshape(xy_data, [], 1);

        % scan_type: 0=raster, 1=CS
        % write it to a .json file
        
        data_name = sprintf('raster_scan_%dpix_%dmic_%.1fHz_fourier%s.json',...
            npix, image_side, raster_freq, rast_mode_name)
        meta_path = fullfile(data_root, data_name)
        
                
        dat= struct('raster_freq', raster_freq,...
            'npix', npix,...
            'width', image_side,...
            'points_per_line', int64(rast.points_per_line),...
            'points_per_period', int64(rast.points_per_line*2),...
            'total_num_points', int64(npix*rast.points_per_line*2),...
            'number_of_scans', 1,...
            'scan_type', 0,...
            'mu_length', image_side*2*npix,...
            'overscan_samples', 0,...
            'pre_scan_samples', 0,...
            'tip_velocity', 0,...
            'actual_sub_samble_perc', 100,...
            'fpga_input', xy_data(:)',...
            'pix_idx', int64([0,0])...
            );
        
        save_local(dat, meta_path);
    end
end

function save_local(data, fpath)
    dest_dir = fileparts(fpath);
    if ~exist(dest_dir, 'file')
        mkdir(dest_dir)
    end
    opts.FileName = fpath;
    savejson('', data,  opts);
end

function y_traj = subl_y_TG(y_height, N)
   if mod(N, 2) ~= 0
       error('Number of samples must be divisible by 2')
   end
   t = (0:N-1)'*AFM.Ts;
   n_by_2 = N/2;
   y_flat = zeros(n_by_2, 1);
   y_up = linspace(0, y_height, n_by_2)';
   y = [y_flat; y_up];
   y_traj = timeseries(y, t);
end



function [freq_even, To_even] = coerce_raster_freq(freq_hz)
    
    To = 1/freq_hz;
    To_even = floor(To/AFM.Ts)*AFM.Ts;
    
    freq_even = 1/To_even;
end

function [u_ff, u_des] = fourier_tri(uu, wo, n_harmonics, Hyr, w_max)
  N = length(uu);
  t = double(0:1:N-1)'*(AFM.Ts);
    
  u_des = 0;
  u_ff = 0;
  for k=1:2:2*n_harmonics
    wk = k*wo;
    if wk > w_max
        break
    end
    fprintf('freq harmonic=%f Hz\n', wk/2/pi);
    phi_k = exp(-1i * wk *t) ; 
    
    ao_des = phi_k' * uu * (2/N); 
    ao_ff = ao_des  /abs((freqresp(Hyr, wk)));

    u_des = u_des + ao_des * phi_k;
    u_ff = u_ff + ao_ff * phi_k;
  end
  u_des = real(u_des + mean(uu));
  u_ff = real(u_ff + mean(uu));

end


