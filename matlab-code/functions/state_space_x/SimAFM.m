classdef SimAFM

  properties
    PLANT;
    controller;
    Nx;
    Nbar;
    useNbar;
    sys_obs;
    L;
    du_max;
    x0_obs;
    x0;
    sys_obs_fp;
    
    thenoise;
    step_amp; % disturbance
    r;
    w;
    rp;
    wp;
    d;
    ws;
    dp;
    wsp;
    gdrift_inv;
    gdrift;
   
    isfxp;
    nw;
    nf;
    simulink_model;
  end

  methods
    
    function self = SimAFM(PLANT, controller, Nx, sys_obs, L, du_max, isfxp, varargin)
    % SimAFM(PLANT, controller, Nx, sys_obs, L, du_max, isfxp, varargin)  
    % varargin is name-value pairs of 
    % 'thenoise', 'step_amp', 'r', 'w', 'rp', 'wp', 'gdrift',
    % 'gdrift_inv', 'nw', 'nf'.
    
      fi_zero = fi(0, 1, 16, 11);
      kd_sys =  tf(1,1, PLANT.Ts);
      
      p = inputParser();
      p.addParameter('thenoise', []);
      p.addParameter('step_amp', 0);
      p.addParameter('r', 0); % These must not be empty, else simulink fails.
      p.addParameter('w', 0);
      p.addParameter('rp', 0);
      p.addParameter('wp', 0);
      p.addParameter('d', 0);
      p.addParameter('ws', 0);
      p.addParameter('dp', fi_zero);
      p.addParameter('wsp', fi_zero);
      p.addParameter('gdrift',  kd_sys);
      p.addParameter('gdrift_inv', kd_sys);
      p.addParameter('nw', 0);
      p.addParameter('nf', 0);
      p.addParameter('useNbar', false)
      p.parse(varargin{:});
      
      self.thenoise = p.Results.thenoise;
      self.step_amp =  p.Results.step_amp;
      self.r = p.Results.r;
      self.w = p.Results.w;
      self.rp = p.Results.rp;
      self.wp = p.Results.wp;
      self.d = p.Results.d;
      self.ws = p.Results.ws;
      self.dp = p.Results.dp;
      self.wsp = p.Results.wsp;
      self.useNbar = p.Results.useNbar;
      
      
      self.PLANT = PLANT;
      self.controller = controller;
      self.Nx = Nx;
      self.Nbar = controller*Nx;

      self.sys_obs = sys_obs;
      self.L = L;
      self.du_max = du_max;
      self.isfxp = isfxp;
      self.x0_obs = sys_obs.b*0;
      self.x0 = PLANT.b*0;
      
      
      self.gdrift = p.Results.gdrift;
      self.gdrift_inv = p.Results.gdrift_inv;
      self.nw = p.Results.nw;
      self.nf = p.Results.nf;
    end
    
    function [Y, U_full, U_nominal, dU, Xhat, Xerr] = sim(sim_obj, ref_traj, dist_traj)
    % NOTE: rather than self, call it sim obj so it is more
    % readable in simulink.  
    % [Y, U, dU] = sim_MPC_fp(self, ref_f)
      %
      % Inputs
      % ------
      %  self : a structure which must have the following feilds
      %
      %    self.K_lqr;
      %    self.PLANT;
      %    self.trun;
      %    self.mpcProb1;
      %    self.du_max;
      %    self.mpc_on;
      %    self.Nx;
      %
      %  ref_f : (scalar) the setpoint to track
      % Outputs
      % -------
      %  Y         : (timeseries) Plant Output
      %  U_full    : U_nominal after running through gd_inv and HSat_inv
      %  U_nominal : (timeseries) Accumulated deltaU
      %  dU : (timeseries) deltaU(k) Control (ie, output of MPC block)
      %  Xhat : state estimation sequence
      %  Xerr : errer state trajectory
      %
      % ------------------------------------------------------------------- %
      % Pull out all the data stored in self to expose it to
      % simulink. There must be a better way...


      fprintf('-----------------------\n');  
      % Expose the sim struct to simulink.
      K_lqr = sim_obj.controller;
      fprintf('simulating as LINEAR\n');

      if ~exist('dist_traj', 'var')
        dist_traj = ref_traj;
        dist_traj.Data = dist_traj.Data*0;
      end
      PLANT = sim_obj.PLANT;
      trun = ref_traj.Time(end);
      
      x0 = SSTools.getNxNu(PLANT)*0;
      uss_0 = 0;
      
      Ts = PLANT.Ts;
      
      % Expose the sim struct to simulink.
      PLANT = sim_obj.PLANT;
      trun = ref_traj.Time(end);
      Ts = PLANT.Ts;
      x0 = sim_obj.x0;
      
      if isempty(sim_obj.thenoise)
        thenoise = ref_traj;
        thenoise.Data = thenoise.Data*0;
      else
        thenoise = sim_obj.thenoise;
      end
      %thenoise = timeseries(ref_traj.Time*0, ref_traj.Time);
      
      sim_obj.step_amp = 0;

      r = sim_obj.r;
      w = sim_obj.w;
      d = sim_obj.d;
      ws = sim_obj.ws;
      
      rp = sim_obj.rp;
      wp = sim_obj.wp;
      dp = sim_obj.dp;
      wsp = sim_obj.wsp;
      
      % sim_obj.gdrift_inv = tf(1,1, PLANT.Ts);
      % sim_obj.gdrift = tf(1,1, PLANT.Ts);
      
      [ ndist, Ns_obs] = size(sim_obj.sys_obs.c);
      Ident_obs = eye(Ns_obs);
      
      % the last row
      C_ydist = Ident_obs(end-ndist+1:end, :);
      % all rows but the last 1
      Ident_obs = Ident_obs(1:end-ndist, :);
      
      % ------------------------  RUN THE SIM ---------------------------- %
      options = simset('SrcWorkspace','current');
      if ~sim_obj.isfxp
        fprintf('simulating as floating-point\n');
        sim('AFMss_fp_obshas_uk', [], options)
      else
        % Turn off warning about fxp states not getting logged.
        fprintf('simulating as FIXED-point\n');
        nw = sim_obj.nw;
        nf = sim_obj.nf;
        id = 'Simulink:Engine:BlkIgnoringUsedAsDStateFlag';
        warning('off', id);
        
        sim('FXP_AFMss_obshas_uk', [], options)
      end
      % provides Y, U, dU, Xhat, U_nominal      
    
    end % sim
   
    function write_control_data(self, data_path, ref, varargin)
    % write_control_data(self, data_path, ref, varargin)
      
      if isempty(self.sys_obs_fp)
        error(['to write the controller data, set property ' ...
               'sys_obs_fp'])
      end
    
      if ~isnumeric(ref)
        ref_traj_path = varargin{1};
        ref_f_1 = 0;
        fid = fopen(ref_traj_path, 'w+');
        fprintf(fid, '%.12f, ', ref.Data(:)'); % as a row.
        fprintf(fid, '\n');
        fclose(fid);
      end

      Ns = size(self.sys_obs.b,1);
      Ns_mpc = Ns; %*size(self.sys_obs.b,2);
      umax = 0;
      Nhyst = length(self.rp);
      Nsat  = length(self.dp);
      Ndrift = length(pole(self.gdrift_inv));
      if Nhyst > 1
        hyst_vec = [self.wp(:); self.rp(:)];
      else
        Nhyst = 0;
      end
      
      if Nsat > 1
        sat_vec = [self.wsp(:); self.dp(:)];
      else
        Nsat = 0;
      end
      if Ndrift >0
        Ndrift_p1 = Ndrift + 1;
        [a,b,c, d] = ssdata(balreal(self.gdrift_inv));
        ABCD = [a, b; c, d];
        % reshape into a column vector along rows.
        ABCD_vec = reshape(ABCD', [], 1);
      else
        Ndrift_p1 = 0;
      end
      K = self.controller;

      AllMatrix = packMatrixDistEst(self.sys_obs_fp,...
        double(self.L), double(K), []);
      
      % ----------- Open File and write data -----------------------------
      [fid, msg] = fopen(data_path, 'w+');
      fprintf(fid, '%d, %f, %.12f,  %.12f, %d, %d, %d\n',...
        Ns, umax, self.du_max, self.Nbar, Nhyst, Nsat, Ndrift_p1);
      
      % fprintf(fid, '0, 0, 0, 0\n'); % was for MPC parameters
        
      % These get written all as one row. 
      fprint_row(fid, '%.12f', AllMatrix);
      if Nhyst > 1
        fprint_row(fid, '%.12f', hyst_vec);
      else
        fprintf(fid, '\n');
      end
      
      if Nsat > 1
        fprint_row(fid, '%.12f', sat_vec);
      else
        fprintf(fid, '\n');
      end
      
      if Ndrift_p1 > 0
        fprint_row(fid, '%.12f', ABCD_vec);
      end
      fclose(fid);
      

    end % write control data
    
    function write_control_data_json(self, data_path, varargin)
    % write_control_data(self, data_path, ref, varargin)
      
      if isempty(self.sys_obs_fp)
        error(['to write the controller data, set property ' ...
               'sys_obs_fp'])
      end
    

      Ns = size(self.sys_obs.b,1);
      umax = 0;
      Nhyst = length(self.rp);
      Nsat  = length(self.dp);
      Ndrift = length(pole(self.gdrift_inv));
      if Nhyst > 1
        hyst_vec = double([self.wp(:)', self.rp(:)']);
      else
        Nhyst = 0;
        hyst_vec = [];
      end
      
      if Nsat > 1
        sat_vec = double([self.wsp(:)', self.dp(:)']);
      else
        Nsat = 0;
        sat_vec = [];
      end
      if Ndrift >0
        Ndrift_p1 = Ndrift + 1;
        [a,b,c, d] = ssdata(balreal(self.gdrift_inv));
        ABCD = [a, b; c, d];
        % reshape into a row vector along rows.
        % We need rows, otherwise JSON parsing in labview will fail.
        ABCD_vec = reshape(ABCD', [], 1)';
      else
        Ndrift_p1 = 0;
        ABCD_vec = [];
      end
      K = self.controller;

      AllMatrix = packMatrixDistEst(self.sys_obs_fp,...
        double(self.L), double(K), []);
      
      % ----------- Open File and write data -----------------------------
      control_data = struct('Ns', Ns, 'umax', umax, 'du_max', double(self.du_max),...
        'Nbar', double(self.Nbar), 'Nhyst', Nhyst, 'Nsat', Nsat, 'Ndrift_p1', Ndrift_p1,...
        'hyst_vec', hyst_vec, 'sat_vec', sat_vec, 'ABCD_vec', ABCD_vec,...
        'AllMatrix_vec', AllMatrix(:)');
      opt.FloatFormat = '%.12f';
      opt.FileName = data_path;
      savejson('', control_data, opt)
      
    end % write control data
   end %methods
end

function fprint_row(fid, fmt, row)
  fprintf(fid, [fmt, ', '], row(1:end-1));
  fprintf(fid, [fmt, '\n'], row(end));
  
end


% allmatrix = packMatrix(sys_obs, L, K, Nx)
%
% Packs observer system matrices, control gain, and observer gain into a
% single column vector. To be read by labview for experiement.
% This is performs the same functionality but in addition to packing 
% [B;L;K';Arows;]
% it appends xss to the end s.t
% [B;L;K';Arows; xss]. To use this with MPC, provde K as empty or zero.
%
% ASSUMES A IS PROVIDED AS A~ = A-LC
function AllMatrix = packMatrixDistEst(sys_obs, L,K, Nx)
    A = sys_obs.a;
    B = sys_obs.b;
    % C = sys_obs.c;
    if isempty(K)
      K = B*0;
    else
      K_T = K';
    end
    
    Ns = size(B, 1);

    AllMatrix = [B(:); 
                 L(:);
                 K_T(:)]; % filler for K
    
    for i=1:Ns
       AllMatrix = [AllMatrix; A(i,:)']; 
    end
    
    AllMatrix = [AllMatrix;
                 Nx];

end