classdef CsExp < handle
  properties
    x;
    x_positive;
    y;
    y_positive;
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
      tmp = load(cs_paths.parent_meta_path); % Provides CsExpMetaIn
      self.meta_in = tmp.CsExpMetaIn;

      self.Ts = Ts;
      self.feature_height = AFM.nm2volts_z*feature_height_nm;
      self.npix = self.meta_in.npix;
      self.width = self.meta_in.width;
      
      self.channel_map = channel_map;
      
      % Get indices for each state.
      self.met_ind = dat_meas(:, channel_map.met);
      % convert the meta cs-measurment index to -4.
      self.met_ind(self.met_ind > 0) = -4;      
      self.idx_state_s = CsExp.divide_by_state(self.met_ind);
      
      self.x = dat_meas(:, channel_map.x);
      self.x = self.x; %- min(self.x); % move to positive orthant.
      self.y = dat_meas(:, channel_map.y);
      self.y = self.y; %- min(self.y); % move to positive orthant.
      
    
      self.t = (0:length(self.x)-1)'*AFM.Ts;
      if exist('gg', 'var') && isa(gg, 'lti')
        self.uz = lsim(gg, (dat_meas(:, channel_map.uz)), self.t);
      else
        self.uz = dat_meas(:, channel_map.uz);
      end
      self.ze = dat_meas(:, channel_map.ze);
      self.ze = dat_meas(:, channel_map.ze);

      
      
      self.Img_raw = zeros(self.npix, self.npix);
      self.Img_smp1d = zeros(self.npix, self.npix);
      self.Img_bp = zeros(self.npix, self.npix);
      self.pix_mask = zeros(self.npix, self.npix);
      self.Gz = zpk([], [], 1, self.Ts);
      
      state_ticks = self.meta_exp.state_counts;
      self.state_times = state_ticks*self.Ts;
      self.time_total = sum(self.state_times);
    end
     
    [sig_psd, freqs, k] = psd_from_intervals(self, signal, state, ...
                                             starts, ends)
    
    function [CS_idx, start_idx, end_idx] = find_cycle_idx(self, time)
    % find the CS-cycle index corresponding to time. A single cycle is defined
    % as xymove --> tip-down--> tip-settle-->scan-->tip-up.
    % very niave search through everything. Would be faster to bisect.

      for CS_idx=1:length(self.idx_state_s.tup)
        % first state is xy-move, last is tip-up. Thus, it is sufficient to check if
        % time is between the first time of move and last time of tup.
        idx_mov = self.idx_state_s.move{CS_idx};
        idx_tup = self.idx_state_s.tup{CS_idx};
        t_mov = self.t(idx_mov);
        t_up = self.t(idx_tup);
        if t_mov(1) <= time && time <= t_up(end)
          start_idx = idx_mov(1);
          end_idx = idx_tup(end);
          return
        end
      end
      % if we get here, we didnt find it.
      CS_idx = [];
      start_idx = [];
      end_idx = [];
      warning('time not found\n');
    end

    function xy_positive(self)
      % Move x-y data to the positive othant. Need to do this based on
      % MEASUREMENT DATA, because with pre-scan, we may be purposefully outside
      % e.g., to the left of the y-axis. Moving everything to positive orthant
      % will leave us with a strip of "unmeasured" data in the final image.
      xmin_meas = [];
      ymin_meas = []; % n.b. min([ [], 9]) = 9.
      for k=1:length(self.idx_state_s.scan)
        idx_k = self.idx_state_s.scan{k};
        xmin_k = min(self.x(idx_k));
        ymin_k = min(self.y(idx_k));
        xmin_meas = min([xmin_meas, xmin_k]); % brackets necessary
        ymin_meas = min([ymin_meas, ymin_k]); % brackets necessary
      end
      
      self.x_positive = self.x - xmin_meas;
      self.y_positive = self.y - ymin_meas;
    end
    
    function plot_all_cycles(self, ax1, ax2, ax3, ax4)
    % plot_time(self, ax1, ax2)
      indc = {'k',        'r', [0, .75, .75], 'b', [.93 .69 .13], ;
       'xy-move', 'tip down', 'tip settle',  '$\mu$-path scan', 'tip up',};
      
      state_seq = {'move', 'tdown', 'tsettle', 'scan', 'tup'};
      hold(ax1, 'on')
      hold(ax2, 'on')
      title(ax1, 'uz')
      title(ax2, 'z-err')
      grid(ax1, 'on')
      grid(ax2, 'on')
      if exist('ax3', 'var')
        hold(ax3, 'on')
        title(ax3, 'x')
        grid(ax3, 'on')
      end
      if exist('ax4', 'var')
        hold(ax4, 'on')
        title(ax4, 'y')
        grid(ax4, 'on')
      end   
      
      Num_cycles = min([length(self.idx_state_s.move), length(self.idx_state_s.tdown),...
        length(self.idx_state_s.tsettle), length(self.idx_state_s.scan),...
        length(self.idx_state_s.tup)]);
      
      for idx_cs_seq = 1:Num_cycles
        
        for k=1:length(state_seq)
          try
          idx_state = self.idx_state_s.(state_seq{k}){idx_cs_seq};
          catch
            keyboard
          end
          uz_k = self.uz(idx_state);
          ze_k = self.ze(idx_state);
          t_k  = idx_state*self.Ts;
          
          plot(ax1, t_k, uz_k, 'color', indc{1, k});
          plot(ax2, t_k, ze_k, 'color', indc{1, k});
          if exist('ax3', 'var')
            x_k = self.x(idx_state);
            plot(ax3, t_k, x_k, 'color', indc{1, k});
          end
          if exist('ax4', 'var')
            y_k = self.y(idx_state);
            plot(ax4, t_k, y_k, 'color', indc{1, k});
          end          
        end
        
      end
      if exist('ax4', 'var')
        linkaxes([ax1, ax2, ax3, ax4], 'x')
      elseif exist('ax', 'var')
        linkaxes([ax1, ax2, ax3], 'x')
      else
        linkaxes([ax1, ax2], 'x')
      end
      plot(ax2, [self.t(1), self.t(end)], [.05, .05], '--k')
      plot(ax2, [self.t(1), self.t(end)], -[.05, .05], '--k')

    end
    

      
    function print_state_times(self)
      tmove = self.state_times(1);
      tlower = self.state_times(2);
      tsettle = self.state_times(3);
      tscan = self.state_times(4);
      tup = self.state_times(5);
      
      fprintf(['Total times\n--------\n',...
               'move   |  lower  |  settle  | scan   | up \n']);
      fprintf('%.3f | %.3f  | %.3f   | %.3f  | %.3f |\n', tmove, tlower, tsettle, tscan, tup);
    end

    function [x_k, y_k, uz_k, ze_k, t_k] = get_state_cycle_k(self, k, state)
      idx_k = self.idx_state_s.(state){k};
      x_k = self.x(idx_k);
      y_k = self.y(idx_k);
      uz_k = self.uz(idx_k);
      ze_k = self.ze(idx_k);
      t_k = idx_k*self.Ts;
      
    end    
    
    function [x_k, y_k, uz_k, ze_k, t_k] = get_tup_k(self, k)
      [x_k, y_k, uz_k, ze_k, t_k] =self.get_state_cycle_k(k, 'tup');
    end
    function [x_k, y_k, uz_k, ze_k, t_k] = get_settle_k(self, k)
      [x_k, y_k, uz_k, ze_k, t_k] =self.get_state_cycle_k(k, 'tsettle');
    end
      
    function [x_k, y_k, uz_k, ze_k, t_k] = get_down_k(self, k)
      [x_k, y_k, uz_k, ze_k, t_k] =self.get_state_cycle_k(k, 'tdown');
    end
      
    function [x_k, y_k, uz_k, ze_k, t_k] = get_scan_k(self, k, positive_xy)
      if exist('positive_xy', 'var') && positive_xy
        idx_k = self.idx_state_s.scan{k};
        x_k = self.x_positive(idx_k);
        y_k = self.y_positive(idx_k);
        [~, ~, uz_k, ze_k, t_k] =self.get_state_cycle_k(k, 'scan');
      else
        [x_k, y_k, uz_k, ze_k, t_k] =self.get_state_cycle_k(k, 'scan');
      end
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

        scan_start_idx = length(U_down) + length(U_settle) + 1;

        t_k = (0:length(U_orig)-1)'*self.Ts;
        U_z = lsim(self.Gz, U_orig(:)-U_orig(1), t_k);
        %U_z = U_orig;    
        U_scan = U_z(scan_start_idx:end);
    end
    
    function [pix_mask] = process_cs_data(self, verbose, figs)
      % bin all the data into pixels. 
      tend_last = 0;
      microns_per_volt = AFM.volts2mic_xy; 
      pix_per_volt = (self.npix/self.width)*microns_per_volt;
      if isempty(self.x_positive) || isempty(self.y_positive)
        self.xy_positive();
      end
      if nargin <2
        verbose = false;
      end
      if verbose && ~exist('figs', 'var')
          figs{1} = figure;
          figs{2} = figure;
          figs{3} = figure;
      end
      if verbose
        [~, ax1] = parse_fig_ax(figs{1});
        [~, ax2] = parse_fig_ax(figs{2});
        [~, ax3] = parse_fig_ax(figs{3});
      end
      
      % Reset the pix_mask
      self.pix_mask = self.pix_mask*0;
      for k = 1:length(self.idx_state_s.scan)
        % Get the data for the current mu-path.
        [X_raw, Y_raw] = self.get_scan_k(k, true);
        Y_raw = Y_raw*pix_per_volt;
        X_raw = X_raw*pix_per_volt;
     
        [U_scan, U_z, U_orig] = self.dynamic_detrend(k);
        
        if max(U_scan) - min(U_scan) > 0.4 % throw out rediculous data.
          fprintf('skipping cycle %d\n', k)
          continue
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




% %     function self = fit_gdrift_per_cycle(self, ax1, ax2, go, gvib)
% %       indc = {'k',        'r', [0, .75, .75], 'b', [.93 .69 .13], ;
% %        'xy-move', 'tip down', 'tip settle',  '$\mu$-path scan', 'tip up',};
% %       
% %       state_seq = {'move', 'tdown', 'tsettle', 'scan', 'tup'};
% %       hold(ax1, 'on')
% %       hold(ax2, 'on')
% % %       title(ax1, 'uz')
% % %       title(ax2, 'z-err')
% %       grid(ax1, 'on')
% %       grid(ax2, 'on')
% %       [z, p, k] = zpkdata(go, 'v');
% %       theta0 = [z;p;k*2];
% %       np = length(p);
% %       ub = ones(2*np+1);
% %       ub(end) = Inf;
% %     
% %       for idx_cs_seq = 1:length(self.idx_state_s.move);
% %         
% %           idx_down = self.idx_state_s.tdown{idx_cs_seq};
% %           idx_settle = self.idx_state_s.tsettle{idx_cs_seq};
% %           idx_scan = self.idx_state_s.scan{idx_cs_seq};
% %           scan_start_idx = length(idx_down) + length(idx_settle) + 1;
% %           
% %           uz_ = self.uz([idx_down, idx_settle]);
% %           ze_ = self.ze([idx_down, idx_settle]);
% %           t_  = [idx_down, idx_settle;]*self.Ts;
% %           
% %           
% %           uz_whole = self.uz([idx_down, idx_settle, idx_scan]);
% %           ze_whole = self.ze([idx_down, idx_settle, idx_scan]);
% %           t_whole  = [idx_down, idx_settle, idx_scan]*self.Ts;
% %           
% %           u0 = uz_(1);
% %           ze0 = ze_(1);
% %           uz_ = uz_ - u0;
% %           ze_ = ze_ - ze0;
% %           
% %           fun = @(theta) fit_gdrift(theta, gvib, ze_, uz_, t_, np);
% %           
% %           theta = lsqnonlin(fun, theta0, 0.9*ub, ub);
% %           gdrift = zpk(theta(np+1:end-1), theta(1:np), theta(end), self.Ts);
% %           z_fit1 = lsim(gdrift*gvib, uz_whole-u0, t_whole)
% %           z_fit2 = -lsim(gdrift*gvib, uz_whole-u0, t_whole) + u0;
% %           
% %           k = idx_down(1);
% %           plot(ax1, t_whole(idx_down-k+1), z_fit2(idx_down-k+1), 'color', indc{1,2}, 'LineStyle', '--')          
% %           plot(ax1, t_whole(idx_settle-k+1), z_fit2(idx_settle-k+1), 'color', indc{1,3}, 'LineStyle', '-')
% %           plot(ax1, t_whole(idx_scan-k+1), z_fit2(idx_scan-k+1), 'color', indc{1,4}, 'LineStyle', '--')
% %           
% %           plot(ax1, t_whole(idx_down-k+1), uz_whole(idx_down-k+1), 'color', indc{1,2}, 'LineStyle', '-')          
% %           plot(ax1, t_whole(idx_settle-k+1), uz_whole(idx_settle-k+1), 'color', indc{1,3}, 'LineStyle', '-')
% %           plot(ax1, t_whole(idx_scan-k+1), uz_whole(idx_scan-k+1), 'color', indc{1,4}, 'LineStyle', '-')
% % 
% %           plot(ax2, t_whole(idx_down-k+1), z_fit1(idx_down-k+1), 'color', indc{1,2}, 'LineStyle', '--')          
% %           plot(ax2, t_whole(idx_settle-k+1), z_fit1(idx_settle-k+1), 'color', indc{1,3}, 'LineStyle', '-')
% %           plot(ax2, t_whole(idx_scan-k+1), z_fit1(idx_scan-k+1), 'color', indc{1,4}, 'LineStyle', '--')
% % 
% %           plot(ax2, t_whole(idx_down-k+1), ze_whole(idx_down-k+1)-ze0, 'color', indc{1,2}, 'LineStyle', '-')          
% %           plot(ax2, t_whole(idx_settle-k+1), ze_whole(idx_settle-k+1)-ze0, 'color', indc{1,3}, 'LineStyle', '-')
% %           plot(ax2, t_whole(idx_scan-k+1), ze_whole(idx_scan-k+1)-ze0, 'color', indc{1,4}, 'LineStyle', '-')
% %           
% %         keyboard
% %       end
% % 
% % %       linkaxes([ax1, ax2], 'x')
% % %       plot(ax2, [self.t(1), self.t(end)], [.05, .05], '--k')
% % %       plot(ax2, [self.t(1), self.t(end)], -[.05, .05], '--k')
% % 
% %     end      
