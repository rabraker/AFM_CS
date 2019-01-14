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
    gg;
    cs_paths;

    UserData = [];
  end

  methods
    function self = CsExp(cs_paths, varargin)
    % channel_map, Ts, feature_height_nm, gg)
    % Obtain cs_paths (which is a struct) from, e.g.,
    % cs_exp_paths(data_root, data_name)
      default_chan_map = ChannelMap([1:5]);
      optP = inputParser();
      optP.addParameter('reload_raw', false, @(s)islogical(s));
      optP.addParameter('channel_map', default_chan_map);
      optP.addParameter('feature_height', Inf);
      optP.addParameter('gg', []);
      optP.addParameter('Ts', AFM.Ts);
      optP.addParameter('load_full', false, @(s)islogical(s));
      optP.parse(varargin{:});
      opts = optP.Results;
      % Check to see if an already processed mat file exists.
      mat_exists = exist(cs_paths.data_path_mat, 'file');

      if ~opts.reload_raw && mat_exists == 2
        fprintf('loading from mat file...\n')
        self = self.load_mat(cs_paths.data_path_mat, opts.load_full);
        fprintf('done\n')
      else
        fprintf('loading from raw data...\n')
        self = self.load_raw_data(cs_paths, opts);
        fprintf('done\n')
      end

    end

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
        save(self.cs_paths.data_path_mat, '-struct', 'self_struct', '-append');
      else
        save(self.cs_paths.data_path_mat, '-struct', 'self_struct');
      end
      warning('on', 'MATLAB:structOnObject');
    end    
    
    function self = load_mat(self, data_path_mat, load_full)
    % Load ourself from the location contained in data_path_mat. If
    % load_full=false (the default), then the original time-series
    % data, x,y,uz,ze, x_positive,y_positive will not be loaded.
    % This is to speed things up when we just want to work with the
    % already processed images.

    % I tried use the matfile. This saves us zero time. I cant
    % figure out how to do this without creating an extra structure
    % and passing that into or out of self. This is wastful...
      self_mat = matfile(data_path_mat);
      to_load_list = properties(self);
      if ~load_full
        no_loads = {'x', 'y', 'x_positive', 'y_positive', 'uz', 'ze'};
        to_load_list = setdiff(to_load_list, no_loads);
      end

      data = load(data_path_mat, to_load_list{:});

      for prop = to_load_list'
        prop = prop{1};
        self.(prop) = data.(prop);
      end

    end

    % Defined in psd_from_intervals.m
    [sig_psd, freqs, k] = psd_from_intervals(self, signal, state, ...
                                             starts, ends)
    % Defined in process_cs_data.m
    [pix_mask] = process_cs_data(self, verbose, figs)
      
    function [CS_idx, start_idx, end_idx] = find_cycle_idx(self, time)
    % find the CS-cycle index corresponding to time. A single cycle is defined
    % as:
    % xymove --> tip-down--> tip-settle-->scan-->tip-up.
    %
    % Very niave search through everything. Would be faster to bisect.

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

    function solve_smp1d(self, recalc)
      if nargin <2
        recalc = false;
      end
      if sum(self.Img_smp1d(:)) ~= 0 && ~recalc
        fprintf('SMP1D already solved. Skipping\n')
        return
      end

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


    function solve_bp(self, recalc, use_2d)
    % Solve the Basis Pursuit problem in either 1d or 2d. If in 1D, use the mex
    % function. 
    % Options
    % -------
    % recalc : (true|false), default false. Do not optimize if self.Img_bp is
    %          non-empty and non-zero.
    % use_2d : (true|false), default false. If true, compute using 2D-dct.
      if nargin <2
        recalc = false;
      end
      if nargin <3
        use_2d = false;
      end
      if ~recalc && ~isempty(self.Img_bp) && sum(self.Img_bp(:)) ~= 0
        warning(['BP solution already calculated, so skipping optimization.',...
          'Pass recalc flag to recompute']);
        return;
      end
      
      [n m] = size(self.Img_raw);
      
      tic
      
      pix_idx = find(CsTools.pixmat2vec(self.pix_mask) > 0.5);
      b = CsTools.pixmat2vec(self.Img_raw);
      b = b(pix_idx);
      
      % y, set of measurements. have to remove all the spots we didn't sample.
      opts = CsTools.l1qc_opts();
      if use_2d
        A = @(x) CsTools.Afun_dct2(x, pix_idx, n);
        At = @(x) CsTools.Atfun_dct2(x, pix_idx, n);
        x0 = At(b);
        eta_vec = CsTools.l1qc_logbarrier(x0, A, At, b, opts);
        self.Img_bp = idct2(CsTools.pixvec2mat(eta_vec, n));
      else
        % A = @(x) CsTools.Afun_dct(x, pix_idx);
        At = @(b) CsTools.Atfun_dct(b, pix_idx, n*m);
        x0 = At(b);
        eta_vec = CsTools.l1qc(x0, b, pix_idx-1, opts);
        self.Img_bp = CsTools.pixvec2mat(idct(eta_vec), n);
      end
      
      time_bp = toc;
      
      fprintf('BP Time: %f\n', time_bp);
      
    end

  end
  
  methods (Access = 'private')
    function flag = raw_data_loaded(self)
      flag = ~isempty(self.x) || ~isempty(self.y)...
        || ~isempty(self.uz) || ~isempty(self.ze);
    end
  end
  
  methods (Static)

    function [figs, axs] = make_traj_figs(figbase)
      Fig_uz = figure(20+figbase); clf

      ax1 = gca();
      Fig_ze = figure(30+figbase); clf
      ax2 = gca();

      Fig_x = figure(40+figbase); clf
      ax3 = gca();
      Fig_y = figure(50+figbase); clf
      ax4 = gca();

      figs = {Fig_uz, Fig_ze, Fig_x, Fig_y};
      axs = {ax1, ax2, ax3, ax4};
    end

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
