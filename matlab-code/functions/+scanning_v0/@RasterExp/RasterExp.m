classdef RasterExp < matlab.mixin.Copyable
% rast_exp = RasterExp(raster_paths, varargin)
% 
% 'rel
  properties
    raster_paths;
    channel_map;
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
    gg;
    
    uz;
    ze;
    x;
    y;
    pix_mat;
    pix_mask;
    UserData;
  end
  
  methods
    function self= RasterExp(raster_paths, varargin)
    % self= RasterExp(raster_paths, varargin)
    
    % [datmat, samps_period, samps_line] 
      self.raster_paths = raster_paths;
      default_chan_map = ChannelMap([1:4, NaN]);
      optP = inputParser();
      optP.addParameter('reload_raw', false, @(s)islogical(s));
      optP.addParameter('channel_map', default_chan_map);
      optP.addParameter('gg', []);
      optP.addParameter('load_full', false, @(s)islogical(s));
      
      optP.parse(varargin{:});
      opts = optP.Results;
      
      % Check to see if an already processed mat file exists.
      mat_exists = exist(raster_paths.data_path_mat, 'file');
      if ~opts.reload_raw && mat_exists
        fprintf('loading from mat file...\n')
        self = self.load_mat(raster_paths.data_path_mat, opts.load_full);
        fprintf('done\n')
      else
        fprintf('loading from raw data...\n')
        self = self.load_raw_data(raster_paths, opts);
        fprintf('done\n')
      end      

    end
    
    % Methods defined in other files
    self = load_raw_data(self, raster_paths, npix, width, opts)
    [ self] = bin_raster_really_slow(self, line_detrender)
    trace_inds = get_trace_inds(self)
    
    function save(self, force_save)
    % Serialize to a .mat file to the location contained
    % in raster_paths.data_path_mat.
    
      if nargin <2
        force_save = false;
      end

      % Remove empty fields, so we don't overwrite data we potentially didn't
      % load with empty.
      if ~force_save && ~self.raw_data_loaded()
        warning(['Not saving data because the raw data is not loaded',...
          'and force_save flag is false. Saving as an append operation',...
          'is very time consuming so is disabled by default.'])
        return
      end
      
      % Go ahead and save it.
      warning('off', 'MATLAB:structOnObject');
      self_struct = struct(self);
      if ~self.raw_data_loaded() % weve been instructed to save anyway
        for fld=fieldnames(self_struct)'
          
          if isempty(self_struct.(fld{1}))
            self_struct = rmfield(self_struct, fld{1});
          end
        end
        save(self.raster_paths.data_path_mat, '-struct', 'self_struct', '-append');
      else
        save(self.raster_paths.data_path_mat, '-struct', 'self_struct');
      end
      warning('on', 'MATLAB:structOnObject');
    end

    function self = load_mat(self, data_path_mat, load_full)
    % Load ourself from the location contained in data_path_mat. If
    % load_full=false (the default), then the original time-series
    % data, x,y,uz,ze, x_positive,y_positive will not be loaded.
    % This is to speed things up when we just want to work with the
    % already processed images.

    % I tried use the matfile() function. This saves us zero time. I cant
    % figure out how to do this without creating an extra structure
    % and passing that into or out of self. This is wastful...
    %%%       self_mat = matfile(data_path_mat);
    
      to_load_list = properties(self);
      if ~load_full
        no_loads = {'x', 'y', 'uz', 'ze'};
        to_load_list = setdiff(to_load_list, no_loads);
      end

      data = load(data_path_mat, to_load_list{:});

      for prop = to_load_list'
        self.(prop{1}) = data.(prop{1});
      end

    end
    
    
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
  
  methods (Access = 'private')
    function flag = raw_data_loaded(self)
      flag = ~isempty(self.x) || ~isempty(self.y)...
        || ~isempty(self.uz) || ~isempty(self.ze);
    end
  end
    
end
