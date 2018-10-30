classdef CsExp < handle
  properties
    x;
    y;
    uz;  
    ze;
    t;
    Ts;
    met_ind;
    idx_state_s;
    npix;
    width;
    channel_map;
    Img_raw;
    Img_smp1d;
    Img_bp;
    pix_mask;
    Gz;
    meta_exp;
    meta_in;
    state_times;
    time_total;
    feature_height
  end
  
  methods
    function self = CsExp(cs_paths, channel_map, Ts, feature_height_nm, gg)
    % Obtain cs_paths (which is a struct) from, e.g.,
    % cs_exp_paths(data_root, data_name)
      
      if ~isa(channel_map, 'ChannelMap')
        error(['channel_map must be of class ChannelMap, but is a ' ...
               '%s'], class(channel_map));
      end
      
      dat_meas = csvread(cs_paths.data_path);
      tmp = load(cs_paths.meta_path);  % Provides ExpMetaData
      self.meta_exp = tmp.ExpMetaData;
      tmp = load(cs_paths.meta_in_path); % Provides CsExpMetaIn
      self.meta_in = tmp.CsExpMetaIn;

      self.Ts = Ts;
      self.feature_height = AFM.nm2volts_z*feature_height_nm;
      self.npix = self.meta_in.npix;
      self.width = self.meta_in.width;
      
      self.channel_map = channel_map;
      self.x = dat_meas(:, channel_map.x);
      self.x = self.x - min(self.x); % move to positive orthant.
      self.y = dat_meas(:, channel_map.y);
      self.y = self.y - min(self.y); % move to positive orthant.
      
      self.t = (0:length(self.x)-1)'*gg.Ts;
      self.uz = lsim(gg, dat_meas(:, channel_map.uz), self.t);
      self.ze = dat_meas(:, channel_map.ze);
      self.met_ind = dat_meas(:, channel_map.met);
      % convert the meta cs-measurment index to -4.
      self.met_ind(self.met_ind > 0) = -4;
      
      % Get indices for each state.
      self.idx_state_s = CsExp.divide_by_state(self.met_ind);
      
      self.Img_raw = zeros(self.npix, self.npix);
      self.Img_smp1d = zeros(self.npix, self.npix);
      self.Img_bp = zeros(self.npix, self.npix);
      self.pix_mask = zeros(self.npix, self.npix);
      self.Gz = zpk([], [], 1, self.Ts);
      
      state_ticks = self.meta_exp.state_counts;
      self.state_times = state_ticks*self.Ts;
      self.time_total = sum(self.state_times);
    end

    function plot_time(self, ax1, ax2)
    % plot_time(self, ax1, ax2)
      indc = {'k',        'r',       [0, .75, .75], 'm', [.93 .69 .13], 'b';
       'xy-move', 'tip down', 'tip settle',  'na', 'tip up', '$\mu$-path scan'};
      plotbyindex(ax1, self.t, self.uz, self.met_ind, indc);
      title(ax1, 'uz')

      hp = plotbyindex(ax2, self.t, self.ze, self.met_ind, indc);
      title(ax2, 'z-err')
      hold on
      plot([self.t(1), self.t(end)], [.05, .05], '--k')
      plot([self.t(1), self.t(end)], -[.05, .05], '--k')
      linkaxes([ax1, ax2], 'x')
    end
    function [x_k, y_k, uz_k, ze_k] = get_settle_k(self, k)
      idx_k = self.idx_state_s.tsettle{k};

      x_k = self.x(idx_k);
      y_k = self.y(idx_k);
      uz_k = self.uz(idx_k);
      ze_k = self.ze(idx_k);
          
    end
      
    function [x_k, y_k, uz_k, ze_k] = get_down_k(self, k)
      idx_k = self.idx_state_s.tdown{k};
        
      x_k = self.x(idx_k);
      y_k = self.y(idx_k);
      uz_k = self.uz(idx_k);
      ze_k = self.ze(idx_k);
          
    end
      
    function [x_k, y_k, uz_k, ze_k] = get_scan_k(self, k)
      idx_k = self.idx_state_s.scan{k};
        
      x_k = self.x(idx_k);
      y_k = self.y(idx_k);
      uz_k = self.uz(idx_k);
      ze_k = self.ze(idx_k);
      
    end
    
    function [U_scan, U_z, U_orig] = dynamic_detrend(self, idx)
        [~, ~, U_scan] = self.get_scan_k(idx);
        [~, ~, U_down] = self.get_down_k(idx);
        [~, ~, U_settle] = self.get_settle_k(idx);
        U_orig = [U_down(:); U_settle(:); U_scan(:)];
        scan_start_idx = length(U_down) + length(U_settle) + 1;
        t_k = (0:length(U_orig)-1)'*self.Ts;
        U_z = lsim(self.Gz, U_orig(:)-U_orig(1), t_k);
        
        U_scan = U_z(scan_start_idx:end);
    end
    function [U_scan, U_z, U_orig] = dynamic_detrend_ze(self, idx)
        [~, ~, U_scan] = self.get_scan_k(idx);
        [~, ~, U_down] = self.get_down_k(idx);
        [~, ~, U_settle] = self.get_settle_k(idx);
        U_orig = [U_down(:); U_settle(:); U_scan(:)];
        U_orig = [U_scan(:)];
%         scan_start_idx = length(U_down) + length(U_settle) + 1;
        scan_start_idx = 1;
        t_k = (0:length(U_orig)-1)'*self.Ts;
        % U_z = lsim(self.Gz, U_orig(:)-U_orig(1), t_k);
        U_z = fft_notch(U_orig(:), 40e-6, 211, 215);
        %U_z = fft_notch(U_z(:), 40e-6, 28, 31);
        %U_z = U_orig;    
        U_scan = U_z(scan_start_idx:end);
    end
    
    function [pix_mask] = process_cs_data(self, verbose, figs)
      % bin all the data into pixels. 
      tend_last = 0;
      microns_per_volt = AFM.volts2mic_xy; 
      pix_per_volt = (self.npix/self.width)*microns_per_volt;
      if verbose
        if ~exist('figs', 'var')
          error('verbose=true but did not recieve figs')
        end
        
        [Fig1, ax1] = parse_fig_ax(figs{1});
        [Fig2, ax2] = parse_fig_ax(figs{2});
        [Fig3, ax3] = parse_fig_ax(figs{3});
      end
      
      for k = 1:length(self.idx_state_s.scan)
        % Get the data for the current mu-path.
        [X_raw, Y_raw] = self.get_scan_k(k);
        Y_raw = Y_raw*pix_per_volt;
        X_raw = X_raw*pix_per_volt;
     
        [U_scan, U_z, U_orig] = self.dynamic_detrend(k);
     
        if max(U_scan) - min(U_scan) > 0.4 % throw out rediculous data.
          %fprintf('skipping\n')
%           continue
        end
        [y_idx, x_idx, U_k] = self.mu_data2pix(X_raw, Y_raw, U_scan);
        self.Img_raw(y_idx, x_idx) = U_k;
        self.pix_mask(y_idx, x_idx) = 1;
        
        if verbose
          t_k = (0:length(U_z)-1)'*self.Ts;
          idx = length(U_z) - length(U_scan);
          % figure(Fig1); 
          plot(ax1, t_k(1:idx)+tend_last, U_orig(1:idx)-U_orig(end), '--g');
          h_f1_uog = plot(ax1, t_k(idx:end)+tend_last, U_orig(idx:end)-U_orig(end), 'k');
          
          plot(ax1, t_k(1:idx)+tend_last, U_z(1:idx)-U_z(end), '--b');
          h_f1_uz = plot(ax1, t_k(idx:end)+tend_last, U_z(idx:end)-U_z(end), 'r');
          tend_last = t_k(end) + tend_last;
          U_ = U_orig(idx:end);
          t_ = (0:length(U_)-1)'*self.Ts;
          % figure(Fig2); 
          h_f2_uog = plot(ax2, t_, U_ - max(U_), 'k');
          h_f2_uz = plot(ax2, t_, U_z(idx:end) - max(U_z(idx:end)), 'r');

          % -------------------
          % ---- visualize ------
          if abs(max(U_k) - min(U_k))> self.feature_height
            % Then we have an edge. 
            cs = 'r';
          else
            cs = 'b';
          end
          plot(ax3, U_k, 'color', cs)
        end
        drawnow();
      end % main loop
      
      if verbose % draw legends
         h_f1_uog.DisplayName = 'original';
         h_f1_uz.DisplayName = 'Dynamic Detrend';
         h_f2_uog.DisplayName = 'original';
         h_f2_uz.DisplayName = 'Dynamic Detrend';
         legend(ax1, [h_f1_uog, h_f1_uz])
         legend(ax2, [h_f2_uog, h_f2_uz])
      end
      
    end % process_cs_data()

    function [y_idx, x_idx, U_k] = mu_data2pix(self, X_raw, Y_raw, U_ks)
      % Make the assumption that the y-data for each path is constant enough.
      % Since we start at the (0,0) corner the xplane, we'll take the floor,
      % and add 1 for 1-based indexing.
      y_idx = min(self.npix-1, floor(mean(Y_raw))) + 1;

      x_spread = max(X_raw) - min(X_raw);
      xpix_start = floor(min(X_raw));
      npix_path_k = ceil(x_spread); % number of bins for this path.
      
      % Now, define a set of bins for the x-direction. Each bin will have a
      % different number of data points. 
      % If we include 0, then there are three points for two pixel bins.
      xbins = linspace(min(X_raw), max(X_raw), npix_path_k+1); 
      U_k = [];
      x_idx = [];
      for jj = 1:npix_path_k
        if xpix_start+jj > self.npix
          continue
        end
        % Get the indeces correspondiing to the current x-data bin.
        ind_x = find(X_raw >= xbins(jj) & X_raw < xbins(jj+1));
        if ~isempty(ind_x) % Avoid errors if it is empty.
                           % Slice out the corresponding height data, and average it
                           % together since this is a single pixel.
          u_pix_jj = mean(U_ks(ind_x)); 
          
          % Collect the uz data into a string of pixels
          U_k = [U_k; u_pix_jj];
          x_idx = [x_idx, xpix_start+jj];
        end
      end
      % Register the control data to zero. We can do this because we are
      % scanning long enough that we are always guaranteed to exit a hole.
      U_k = U_k - max(U_k);
      
    end % bin_data_into_pix
    
    function solve_smp1d(self)
      [n,m] = size(self.Img_raw);

      tic
      samp_frac = sum(sum(self.pix_mask))/(n*m);
      reduce_frac = 0.1;  %It's my understanding Yufan says 1/10 of sub-sample fraction.
      maxiter = round(samp_frac*reduce_frac*n*m);
      fprintf('sample fraction: %.2f     max-iter: %d\n', samp_frac, maxiter);
      self.Img_smp1d = SMP_1D(self.Img_raw, self.pix_mask, maxiter);
      time_smp = toc;      
      fprintf('Total smp-1d solution time: %.1f\n', time_smp);
    end
    
    function solve_basis_pursuit(self)
      [n m] = size(self.Img_raw);
      tic
      I_vector = PixelMatrixToVector(self.Img_raw);

      pix_mask_vec = PixelMatrixToVector(self.pix_mask);
      I_vector = I_vector(find(pix_mask_vec>0.5));

      A = @(x) IDCTfun(x,pix_mask_vec);
      At = @(x) DCTfun(x,pix_mask_vec);

      Ir_bp = idct(l1qc_logbarrier(At(I_vector), A, At, I_vector, 0.1));
      Ir_bp = real(Ir_bp);
      self.Img_bp = PixelVectorToMatrix(Ir_bp,[n m]);
      time_bp = toc;
      fprintf('BP Time: %f\n', time_bp);

      self.Img_bp = detrend_plane(self.Img_bp);
      
      
    end
    
    
  end
  
  methods (Static)
    
    function idx_state_s = divide_by_state(met_ind)
      names = {'move', 'tdown', 'tsettle', 'scan', 'tup'};
      names_idx = [1:5];
      idx_state_s.scan = {};
      idx_state_s.move = {};
      idx_state_s.tdown = {};
      idx_state_s.tsettle = {};
      idx_state_s.tup = {};      
      
      met_ind_temp = abs(met_ind);
      idx_end = 0;
      
      break_out = false;
      k=1;
      while true 
        
        for j = 1:length(names)
          idx_end_prev = idx_end;
          j_next = max(mod(j+1, 6), 1); % wrap back to 1 when j=5.
          idx_start = find(met_ind_temp(idx_end_prev+1:end) == names_idx(j), 1, 'first') + idx_end_prev;
          
          idx_end   = find(met_ind_temp(idx_end_prev+1:end) == names_idx(j_next), 1, 'first') + idx_end_prev-1;
          
          if isempty(idx_end)
            idx_end   = find(met_ind_temp(idx_end_prev+1:end) == names_idx(j), 1, 'last') + idx_end_prev;
          end
          
          if isempty(idx_end) || idx_end >=length(met_ind_temp)
            break_out = true;
            break
          end
          idx_state_s.(names{j}){k} = idx_start:idx_end;
        end
        
        if break_out
          break
        end
        k = k+1;
      end      
    end
    
     
        
    
    
  end
  
  
end
