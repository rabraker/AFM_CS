

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
    
    volts_per_sample;
%     mu_lengths_mic;
%     mu_lengths_pix;
%     mu_lengths_volts;
    
    mu_length_mic_nom;
    mu_length_volts;
    mu_pix_nom;
    
    XR_pix_starts;
    YR_pix_starts;
    XR_volt_starts;
    YR_volt_starts;
    
    mu_path_traj_s;
    
    pre_pad_samples;
    overscan_samples;
    
    sub_sample_perc;
    met_idx_pad = 10;
    is_connected=false;
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
      p.addParameter('ax', []);
      p.parse(varargin{:});
      self.overscan_samples = p.Results.overscan_samples;
      self.pre_pad_samples =  p.Results.pre_pad_samples;
      ax = p.Results.ax;
      
      self.pix_mask = pix_mask;
      self.npix = size(pix_mask, 2);
      self.width_mic = width_mic;
      self.pix2mic = self.width_mic/self.npix;
      self.pix2volts = (self.width_mic/self.npix) * AFM.mic2volt_xy;
      
      self.x_rate_mic_per_sec = x_rate;
      self.x_rate_pix_per_sec = (self.npix/self.width_mic)*self.x_rate_mic_per_sec;
      self.mu_length_mic_nom = mu_len;
      self.mu_pix_nom = ceil(mu_len * (self.npix/self.width_mic));
      self.mu_length_volts = self.mu_pix_nom * self.pix2volts;
      
      self.Ts = Ts;
      self.volts_per_sample = self.x_rate_mic_per_sec * AFM.mic2volt_xy * self.Ts;
      % self.build_xr_yr_volt_starts();
      self.build_mu_path_traj_s(ax)
      N_paths = length(self.XR_pix_starts);
      
      self.sub_sample_perc = 100 * sum(self.pix_mask(:)) / self.npix^2;
    end
    function fname = get_fname(self)
      perc = floor(self.sub_sample_perc);
      prescan_str = sprintf('%dprescan', self.pre_pad_samples);
      if self.is_connected
        con_str= 'isconnect';
      else
        con_str = 'notconnect';
      end
      if self.N_paths <= 2
        con_str = 'istotallyconnected';
      end
      fname = sprintf('cs-traj-%dpix-%dperc', self.npix, perc);
      fname = sprintf('%s-%dnm-%dmic', fname, self.mu_length_mic_nom*1000, self.width_mic);
      fname = sprintf('%s-%.2dHz', fname, self.x_rate_mic_per_sec / (self.width_mic*2));
      fname = sprintf('%s-%s-%s.json', fname, prescan_str, con_str);
    end

    function write_data_json(self, json_fname)
      vec = self.as_vector();
      % For full compatibility with doing raster as CS, we need extra fields.
      
      % scan_type: 0=raster, 1=CS
      if self.is_connected
        N_paths_effective = self.N_paths + self.met_idx_pad - 1;
      else
        N_paths_effective = self.N_paths + 1 - 1;
      end

      data = struct(...
              'raster_freq', 0,...
              'npix', self.npix,...
              'width', self.width_mic,...
              'points_per_line', int64(0),...
              'points_per_period', int64(0),...
              'total_num_points', int64(length(vec)/3),... % /3 because have x,y,idx
              'number_of_scans', N_paths_effective,...
              'scan_type', 1,...
              'mu_length', self.mu_length_mic_nom,...
              'overscan_samples', self.overscan_samples,...
              'pre_scan_samples', self.pre_pad_samples,...
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
      
      if isempty(self.mu_path_traj_s)
        self.build_mu_path_traj_s();
      end
      
      vec = [];
      % The easiest way to interleve the data suitible for pumping into the
      % FPGA is to build a matrix and then reshape.
      for k=1:length(self.mu_path_traj_s)
        mptc_k = self.mu_path_traj_s{k};
        xt = mptc_k.xt;
        yt = mptc_k.yt;
        met_idx = mptc_k.met_idx;
        try
          if self.is_connected
            met_idx_pad = self.met_idx_pad;
          else
            met_idx_pad = 1;
          end
        assert(all(abs(met_idx) == k + met_idx_pad - 1)); % should only have k and -k.
        catch
          keyboard
        end
        [xrp, yrp] = self.adjust_pre_pad(xt(1), yt(1), length(xt));
        
        xr_pre = linspace(xrp, xt(1), self.pre_pad_samples);
        yr_pre = xr_pre*0 + yrp;
        met_idx_pre = xr_pre*0 + self.met_idx_pad - 1;
        
        % the setpoint has a meta-idx=0;
        vec_k = [[xr_pre(1), xr_pre(:)'];
                 [yr_pre(1), yr_pre(:)'];
                 [0, met_idx_pre(:)']];

        
        met_idx(end) = -1;
        vec_k = [vec_k, [xt(:)'; yt(:)'; met_idx(:)']];  %#ok<AGROW>
% debug the bizzar trajectories we are getting.
%         figure(200)
%         plot(vec_k(1, :))
%         title('x')
%         figure(201)
%         plot(vec_k(2, :))
%         title('y')
%         keyboard        
%                 
%         keyboard        
        
        vec = [vec; reshape(vec_k, [], 1)];              %#ok<AGROW>
      end
      
    end
    
    function mu_path_connect_rad(self, mptc_opts)
      if isempty(self.mu_path_traj_s)
        self.build_mu_path_traj_s();
      end

      mptc_c = {};
      mptc_s = self.mu_path_traj_s;
      for j=1:mptc_opts.Npasses
        k = 1;
        while ~isempty(mptc_s)
          mptc = mptc_s{1};
          
          mptc_s(1) = [];
          
          [mptc_s, mptc] = mpt_connect_rad(mptc_s, mptc, mptc_opts);
          % corerce the meta index to the actual k.
          mptc.met_idx(mptc.met_idx >0) = k + self.met_idx_pad - 1;
          mptc.met_idx(mptc.met_idx <0) = -k - self.met_idx_pad + 1;
          
          mptc_c{k} = mptc;
          
          k=k+1;
        end
      end
      
      self.mu_path_traj_s = mptc_c;
      self.is_connected = true;
      self.N_paths = length(self.mu_path_traj_s);
    end
    
    function [xr, yr, N] = adjust_pre_pad(self, xr_volts, yr_volts, N_)
    % [xr, yr, N] = adjust_pre_pad(self, xr_volts, yr_volts, N_)
    % 
      N = self.pre_pad_samples; %+N_
      
      volts_per_sec = self.x_rate_mic_per_sec * AFM.mic2volt_xy();
      
      x_diff_ticks = N * volts_per_sec * self.Ts; % Ts is seconds per tick
      
      xr = xr_volts - x_diff_ticks;
      yr = yr_volts;
      
    end

    function build_mu_path_traj_s(self,ax)
    % build_mu_path_traj_s(self, verbose)
      if ~exist('ax', 'var')
        ax=[];
      end
      
      N_per_mupath = ceil(self.mu_length_volts/self.volts_per_sample) + self.overscan_samples;
      
      scan_count = 0;
      for j=1:size(self.pix_mask,1)
        xr_idx = self.get_mu_starts(self.pix_mask(j,:));
        
        for k=1:length(xr_idx)
          scan_count = scan_count + 1;
          x1s = xr_idx(k) * self.pix2volts;
          y1s = j * self.pix2volts;
          
          xt = x1s + (0:N_per_mupath-1)' * self.volts_per_sample;
          yt = ones(N_per_mupath, 1) * y1s;
          
          met_idx = ones(N_per_mupath,1)*scan_count;
          
          self.mu_path_traj_s{scan_count} = MuPathC(xt, yt, met_idx);
          
          if ~isempty(ax)
            plot(ax, self.mu_path_traj_s{scan_count}.xt,...
              self.mu_path_traj_s{scan_count}.yt, '-b', 'LineWidth', 1.5);
          end
        end
      end
      self.N_paths = length(self.mu_path_traj_s);
    end


%     function build_xr_yr_pix_starts(self)
%     % build_xy_yr_starts(self)
%       xr = [];
%       yr = [];
%       for j=1:size(self.pix_mask,1)
%         xr_idx = self.get_mu_starts(self.pix_mask(j,:));
% 
%         if ~isempty(xr_idx)
%           yr = [yr; repmat(j, length(xr_idx), 1)];
%           xr = [xr; xr_idx(:)];
%         end
%         
%       end
%       % convert from pixel locations to volts
%       self.XR_pix_starts = xr; %*self.pixel2volts;
%       self.YR_pix_starts = yr; %*self.pixel2volts;
%       self.N_paths = length(xr);
%     end
    
%     function build_xr_yr_volt_starts(self)
%       if isempty(self.XR_pix_starts) || isempty(self.YR_pix_starts)
%         self.build_xr_yr_pix_starts();
%       end
%       
%       self.XR_volt_starts = self.XR_pix_starts * self.pix2volts;
%       self.YR_volt_starts = self.YR_pix_starts * self.pix2volts;
%     end
    

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

% THIS IS THE OLD, ROW CONNECTION CODE THAT DOESNT WORK.
% % %    function connect_mu_paths(self, Tmu)
% % %       % connect_mu_paths(self, Tmu)
% % %       % Connect mu-paths in a single row if the distance separating the end of
% % %       % one and the beginning of the next is less than: 
% % %       %          pix_min = Tmu*self.x_rate_pix_per_sec
% % %       %
% % %       % Updates the pix_mask property
% % %       %
% % %       % Tmu is the time.
% % %       
% % %       % minimum number of pixels needed to separate a mupath for them not to
% % %       % get connected.
% % %       pix_min = Tmu*self.x_rate_pix_per_sec;
% % %       self.mu_lengths_pix = []; % reset to empty.
% % %       
% % %       for j=1:size(self.pix_mask,1)
% % %         xr_idx = self.get_mu_starts(self.pix_mask(j,:));
% % %         if isempty(xr_idx)
% % %           continue
% % %         end
% % %         mu_lens_j = repmat(self.mu_pix_nom, length(xr_idx), 1);
% % %         [xr_idx, mu_lens] = self.connect_row(xr_idx, mu_lens_j, pix_min);
% % %         if length(xr_idx) >1
% % %           fprintf('connecting xr=%d to xr=%d\n', xr_idx(1), xr_idx(2));
% % %         end
% % %         
% % %         self.mu_lengths_pix = [self.mu_lengths_pix; mu_lens];
% % %         for k=1:length(xr_idx)
% % %           self.pix_mask(j, xr_idx(k):xr_idx(k)+mu_lens(k)-1)=1;
% % % %           if length(self.pix_mask(j,:)) > 16
% % % %             keyboard
% % % %           end
% % %         end
% % %       end
% % %       self.mu_lengths_mic = self.mu_lengths_pix * self.pix2mic;
% % %       self.mu_lengths_volts = self.mu_lengths_pix * self.pix2volts;
% % %       self.N_paths = length(self.mu_lengths_pix);
% % %     end
% % %     
% % %     
% % %     function [xr_idx, mu_pix_s] = connect_row(self, xr_idx, mu_pix_s, pix_min)
% % %     % [xr_idx, mu_pix_s] = connect_row(self, xr_idx, mu_pix_s, pix_min)
% % %     % This function works recursively.
% % %       if length(xr_idx) < 2
% % %         return;
% % %       end
% % %       
% % %       if xr_idx(1)+mu_pix_s(1)-1 >= xr_idx(2)
% % %         % They are close enough so connect them.
% % %         mu_pix_s(1) = xr_idx(2) - xr_idx(1) + mu_pix_s(2)+1;
% % %         xr_idx(2) = [];
% % %         mu_pix_s(2) = [];
% % %         
% % %         [xr_idx, mu_pix_s] = self.connect_row(xr_idx(1:end), mu_pix_s(1:end), pix_min);
% % % 
% % %       else
% % %         % otherwise, try to connect the rest of the row.
% % %         [xr_idx_remain, mu_pix_remain] = self.connect_row(xr_idx(2:end), mu_pix_s(2:end), pix_min);
% % %         xr_idx = [xr_idx(1); xr_idx_remain];
% % %         mu_pix_s = [mu_pix_s(1); mu_pix_remain];
% % %       end
% % % 
% % %     end








