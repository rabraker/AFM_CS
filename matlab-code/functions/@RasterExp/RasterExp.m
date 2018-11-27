classdef RasterExp <handle
  
  
  properties
    paths;
    % parent_data;
    xref;
    yref;
    meta_data;
    Ts;
    npix;
    width;
    volts2pix;
    micron2pix;

    samps_per_period;
    samps_per_line;
    
    uz;
    ze;
    x;
    y;
    pix_mat;
    pix_mask;
  end
  
  methods
    function self= RasterExp(raster_paths, npix, width)
    % [datmat, samps_period, samps_line] 
      self.paths = raster_paths;
      
      fprintf('Loading Meta file...\n%s\n', raster_paths.meta_path)
      meta_data = load(raster_paths.meta_path);
      self.meta_data = meta_data.Cluster;
      self.Ts = meta_data.Cluster.raster_scan_params.TsTicks/ ...
                (40e6);
      
      fprintf('Loading Data file...\n%s\n', raster_paths.data_path)
      datmat = csvread(raster_paths.data_path);
      
      % Load parent data.
      parent_dat = csvread(raster_paths.parent_path);
      xyref = reshape(parent_dat', 2, [])';
      self.xref = xyref(:,1);
      self.yref = xyref(:,1);
      % Compute unit conversions
      micron2pix = npix/width;
      volts2pix = AFM.volts2mic_xy  * micron2pix;
      
      self.npix = npix;
      self.volts2pix = volts2pix;
      self.micron2pix = micron2pix;
      self.samps_per_period = size(parent_dat,1)/2; % twice as many in here for x & y.
      self.samps_per_line = self.samps_per_period/2;
      if floor(self.samps_per_line)~= self.samps_per_line
        warning('Non integer number of samples per line = %f.', ...
                self.samps_per_line)
      end
      
      % Pull out data. First, drop anything extra that got collected.
      datmat = datmat([1:npix*self.samps_per_period], :);
      self.x = datmat(:, 1);
      self.y = datmat(:, 2);
      self.ze = datmat(:, 3);
      self.uz = datmat(:, 4);
    end
    %%
    [ self] = bin_raster_really_slow(self, line_detrender)
    
    trace_inds = get_trace_inds(self)
    
    function plot_n_periods(self, signal, ax, N0, N1)
    % plot_n_periods(self, signal, ax, N0, N1)
    % 
    % Plot the signal ('x', 'y', 'ze', 'uz') to axis ax between
    % raster periods N0 and N1.
      sig = self.(signal)(N0 * self.samps_per_line + 1: N1* ...
                          self.samps_per_line);
      
      plot(ax, sig);
    end


    function [uz_row, x_row] = get_row(self, row)
    % [uz_row, x_row] = get_row(row)
    % Extract from the raw data a single row of control and
    % x-direction data, which should
    % correspond to an actual row in the post-processed image.      
      start_idx = self.samps_per_period*(row-1) + 1;
      end_idx = start_idx + self.samps_per_line;

      x_row = self.x(start_idx:end_idx);
      uz_row = detrend(self.uz(start_idx:end_idx) );
      x_row = x_row - x_row(1);
      
      x_row = x_row * (self.npix/x_row(end) );
    end
  end
end
