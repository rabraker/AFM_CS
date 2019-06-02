% Generate x-y input waveforms.
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
lines = [512, 128, 64];
freqs = [1, 2.5, 5, 8, 10, 15];
for pp = 1:length(lines)
%     npix = 128;
    npix = lines(pp);
    for k=1:length(freqs)
        % raster_freq =
        raster_freq = coerce_raster_freq(freqs(k)); % Hz.
        image_side = 5; % micro-meters.

        rast = RasterTraj(image_side, raster_freq, npix);
        
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
        xy_data = zeros(2*length(x_traj.Time), 1);
        
        j = 1;
        for k=1:2:length(xy_data)
            xy_data(k) = r_ff(j);
            xy_data(k+1) = y_traj.Data(j);
            j = j+1;
        end
        
        % scan_type: 0=raster, 1=CS
        % write it to a .json file
        
        data_name = sprintf('raster_scan_%dpix_%dmic_%.1fHz_fourier.json',...
            npix, image_side, raster_freq)
        
        target_dir = sprintf('%dmicrons', image_side);
        data_root = PATHS.raster_image_data('5microns', 'parents-loop');
        
        meta_path = fullfile(data_root, data_name)
        
        opts.FileName = meta_path;
        
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
        
        savejson('', dat,  opts);
    end
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


