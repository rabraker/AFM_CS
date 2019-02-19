

classdef MuPathTraj < handle
  
  properties
    width_mic;
    pix_mask;
    x_rate_mic_per_sec;
    x_rate_pix_per_sec;
    npix;
    N_paths;
    Ts;
    pix2mic;
    pix2volts;
    
    mu_lengths_mic;
    mu_lengths_pix;
    mu_lengths_volts;
    
    mu_length_mic_nom;
    mu_pix_nom;
    
    XR_pix_starts;
    YR_pix_starts;
    XR_volt_starts;
    YR_volt_starts;
    
    pre_pad_samples;
    overscan_samples;
    
    sub_sample_perc;
  end
  

  methods
    function self = MuPathTraj(pix_mask, width_mic, mu_len, x_rate, Ts, varargin)
    % self = MuPathTraj(pix_mask, width_mic, mu_len, x_rate, Ts, varargin)
    %
    % varargin:
    % ---------
    %   'pre_pad_samples', (int) number of samples to move the
    %    trajectory left by (so we can begin the scan during z-settling).
    %    
    %   'overscan_samples': (int) number of samples to overscan by (to correct
    %   for ramp following error).
    
      p = inputParser();
      p.addParameter('pre_pad_samples', 0);
      p.addParameter('overscan_samples', 0);
      p.parse(varargin{:});
      self.overscan_samples = p.Results.overscan_samples;
      self.pre_pad_samples =  p.Results.pre_pad_samples;
      
      self.pix_mask = pix_mask;
      self.npix = size(pix_mask, 2);
      self.width_mic = width_mic;
      self.pix2mic = self.width_mic/self.npix;
      self.pix2volts = (self.width_mic/self.npix) * AFM.mic2volt_xy;
      
      self.x_rate_mic_per_sec = x_rate;
      self.x_rate_pix_per_sec = (self.npix/self.width_mic)*self.x_rate_mic_per_sec;
      self.mu_length_mic_nom = mu_len;
      self.mu_pix_nom = ceil(mu_len * (self.npix/self.width_mic));
      self.Ts = Ts;
      
      self.build_xr_yr_volt_starts();
      N_paths = length(self.XR_pix_starts);
      self.mu_lengths_pix = repmat(self.mu_pix_nom, N_paths, 1);
      self.mu_lengths_mic = self.mu_lengths_pix*self.pix2mic;
      self.mu_lengths_volts = self.mu_lengths_pix * self.pix2volts;
      
      self.sub_sample_perc = 100 * sum(self.pix_mask(:)) / self.npix^2;
    end
  
    function write_data_json(self, json_fname)
      vec = self.as_vector();
      % For full compatibility with doing raster as CS, we need extra fields.
      
      % scan_type: 0=raster, 1=CS
      data = struct(...
              'raster_freq', 0,...
              'npix', self.npix,...
              'width', self.width_mic,...
              'points_per_line', int64(0),...
              'points_per_period', int64(0),...
              'total_num_points', int64(length(vec)/3),... % /3 because have x,y,idx
              'number_of_scans', self.N_paths,...
              'scan_type', 1,...
              'mu_length', self.mu_length_mic_nom,...
              'overscan_samples', self.overscan_samples,...
              'pre_pad_samples', self.pre_pad_samples,...
              'tip_velocity', self.x_rate_mic_per_sec,...
              'actual_sub_samble_perc', self.sub_sample_perc,...
              'fpga_input', vec(:)',...
              'pix_idx', int64(find(self.pix_mask > 0.5))'...
              );
    
     opts.FileName = json_fname;
     savejson('', data, opts)

    end
       
    
    function write_data(self, csv_fname)
      vec = self.as_vector();
      
      csvwrite(csv_fname, vec);
      
      meta_in_fname = strrep(csv_fname, '.csv', '.mat');
      CsExpMetaIn.width = self.width_mic;
      
      CsExpMetaIn.mu_length = self.mu_length_mic_nom;
      CsExpMetaIn.tip_velocity = self.x_rate_mic_per_sec;
      CsExpMetaIn.npix = self.npix;
      CsExpMetaIn.pix_mask = self.pix_mask;
      CsExpMetaIn.actual_perc = self.sub_sample_perc;
      
      save(meta_in_fname, 'CsExpMetaIn');
    end
    
    function vec = as_vector(self, pre_pad_samples)
    % vec = as_vector(self,  pre_pad_samples)
    
      if exist('pre_pad_samples', 'var')
        self.pre_pad_samples = pre_pad_samples;
      end
      
      if isempty(self.XR_volt_starts) || isempty(self.YR_volt_starts)
        self.build_xr_yr_volt_starts();
      end
      self.mu_lengths_volts = self.mu_lengths_pix * self.pix2volts;
      x_rate_volts = self.x_rate_mic_per_sec * AFM.mic2volt_xy();
      x_rate_volts_per_sample = x_rate_volts*self.Ts;
      N_vec = ceil(self.mu_lengths_volts./x_rate_volts_per_sample);
      
      vec = [];
      % The easiest way to interleve the data suitible for pumping into the
      % FPGA is to build a matrix and then reshape.
      for k=1:length(self.XR_volt_starts)
        
        xr_ = self.XR_volt_starts(k);
        yr_ = self.YR_volt_starts(k);
        N_ = N_vec(k);
        [xr_k, yr_k, N_k] = self.adjust_pre_pad(xr_, yr_, N_);
        N_k = N_k + self.overscan_samples;
        
        % the setpoint has a meta-idx=0;
        vec_k = [xr_k; 
                 yr_k;
                 0];
        % Value of x at end of ramp.
        x_N = xr_k + (N_k-1) * x_rate_volts_per_sample;
        x_mu_ramp_k =  linspace(xr_k, x_N, N_k);
        y_mu_k = x_mu_ramp_k * 0 + yr_k;
        
        met_idx = ones(1, N_k) * k;
        met_idx(end) = -1;
        
        vec_k = [vec_k, [x_mu_ramp_k; y_mu_k; met_idx]]; %#ok<AGROW>
        
        vec = [vec; reshape(vec_k, [], 1)];              %#ok<AGROW>
      end
      
    end

    function [xr, yr, N] = adjust_pre_pad(self, xr_volts, yr_volts, N_)
    % [xr, yr, N] = adjust_pre_pad(self, xr_volts, yr_volts, N_)
    % 
      N = N_ + self.pre_pad_samples;
      
      volts_per_sec = self.x_rate_mic_per_sec * AFM.mic2volt_xy();
      
      x_diff_ticks = N_ * volts_per_sec * self.Ts; % Ts is seconds per tick
      
      xr = xr_volts - x_diff_ticks;
      yr = yr_volts;
      
    end
    
    
    function connect_mu_paths(self, Tmu)
      % connect_mu_paths(self, Tmu)
      % Connect mu-paths in a single row if the distance separating the end of
      % one and the beginning of the next is less than: 
      %          pix_min = Tmu*self.x_rate_pix_per_sec
      %
      % Updates the pix_mask property
      %
      % Tmu is the time.
      
      % minimum number of pixels needed to separate a mupath for them not to
      % get connected.
      pix_min = Tmu*self.x_rate_pix_per_sec;
      self.mu_lengths_pix = []; % reset to empty.
      
      for j=1:size(self.pix_mask,1)
        xr_idx = self.get_mu_starts(self.pix_mask(j,:));
        if isempty(xr_idx)
          continue
        end
        mu_lens_j = repmat(self.mu_pix_nom, length(xr_idx), 1);
        [xr_idx, mu_lens] = self.connect_row(xr_idx, mu_lens_j, pix_min);
        
        self.mu_lengths_pix = [self.mu_lengths_pix; mu_lens];
        for k=1:length(xr_idx)
          self.pix_mask(j, xr_idx(k):xr_idx(k)+mu_lens(k)-1)=1;
%           if length(self.pix_mask(j,:)) > 16
%             keyboard
%           end
        end
      end
      self.mu_lengths_mic = self.mu_lengths_pix * self.pix2mic;
      self.mu_lengths_volts = self.mu_lengths_pix * self.pix2volts;
      self.N_paths = length(self.mu_lengths_pix);
    end
    
    
    function [xr_idx, mu_pix_s] = connect_row(self, xr_idx, mu_pix_s, pix_min)
    % [xr_idx, mu_pix_s] = connect_row(self, xr_idx, mu_pix_s, pix_min)
    % This function works recursively.
      if length(xr_idx) < 2
        return;
      end
      
      if xr_idx(1)+mu_pix_s(1)-1 >= xr_idx(2)
        % They are close enough so connect them.
        mu_pix_s(1) = xr_idx(2) - xr_idx(1) + mu_pix_s(2)+1;
        xr_idx(2) = [];
        mu_pix_s(2) = [];
        
        [xr_idx, mu_pix_s] = self.connect_row(xr_idx(1:end), mu_pix_s(1:end), pix_min);

      else
        % otherwise, try to connect the rest of the row.
        [xr_idx_remain, mu_pix_remain] = self.connect_row(xr_idx(2:end), mu_pix_s(2:end), pix_min);
        xr_idx = [xr_idx(1); xr_idx_remain];
        mu_pix_s = [mu_pix_s(1); mu_pix_remain];
      end

    end
    
    function build_xr_yr_pix_starts(self)
    % build_xy_yr_starts(self)
      xr = [];
      yr = [];
      for j=1:size(self.pix_mask,1)
        xr_idx = self.get_mu_starts(self.pix_mask(j,:));
        
        if ~isempty(xr_idx)
          yr = [yr; repmat(j, length(xr_idx), 1)];
          xr = [xr; xr_idx(:)];
        end
        
      end
      % convert from pixel locations to volts
      self.XR_pix_starts = xr; %*self.pixel2volts;
      self.YR_pix_starts = yr; %*self.pixel2volts;
      self.N_paths = length(xr);
    end
    
    function build_xr_yr_volt_starts(self)
      if isempty(self.XR_pix_starts) || isempty(self.YR_pix_starts)
        self.build_xr_yr_pix_starts();
      end
      
      self.XR_volt_starts = self.XR_pix_starts * self.pix2volts;
      self.YR_volt_starts = self.YR_pix_starts * self.pix2volts;
    end
    

  end
  
  methods (Access=private)

    function [xr_idx, mu_len_s] = get_mu_starts(self, mask_row)
    % [xr_idx, mu_len_s] = get_mu_starts(self, mask_row)
    
      % Note that diff will jump at the index before, so that
      % diff([0 0 0 1 1 1 0 0]) = [0 0 1 0 0 -1 0]. -1 is the jump down, from being
      % sampled to not. So we only care about the positive diff, ie, ==1
      
      % If the first column has a mu-path, it wont show up when we take the
      % diff.
      if mask_row(1) == 1
        xr_idx = 1;
      else
        xr_idx = [];
      end
      xr_idx = [xr_idx, find(diff(mask_row) == 1) + 1];
    end
    
   end

end











