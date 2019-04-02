

classdef MuPathBounce < handle
  properties
    N_paths;
    Ts=AFM.Ts;
    
    XR_volt_starts;
    YR_volt_starts;
    samps_per_bounce;
    
    pre_pad_samples=0;
  end
  methods
    function self = MuPathBounce(N_paths, samps_per_bounce, varargin)
      p = inputParser();
      p.addParameter('pre_pad_samples', 0);
      p.parse(varargin{:});
      
      self.pre_pad_samples =  p.Results.pre_pad_samples;
      self.N_paths = N_paths;
      self.XR_volt_starts = zeros(N_paths, 1);
      self.YR_volt_starts = zeros(N_paths, 1);
      self.samps_per_bounce = samps_per_bounce;
    end
    function write_data_json(self, json_fname)
      vec = self.as_vector();
      % For full compatibility with doing raster as CS, we need extra fields.
      % scan_type: 0=raster, 1=CS
      data = struct(...
              'raster_freq', 0,...
              'npix', 0,...
              'width', 0,...
              'points_per_line', int64(0),...
              'points_per_period', int64(0),...
              'total_num_points', int64(length(vec)/3),... % /3 because have x,y,idx
              'number_of_scans', self.N_paths,...
              'scan_type', 1,...
              'mu_length', 0,...
              'overscan_samples', 0,...
              'pre_scan_samples', self.pre_pad_samples,...
              'tip_velocity', 0,...
              'actual_sub_samble_perc', 0,...
              'fpga_input', vec(:)',...
              'pix_idx', int64([0,0])...
              );
    
     opts.FileName = json_fname;
     savejson('', data, opts)

    end
    
    function vec = as_vector(self)
      vec = [];
      % The easiest way to interleve the data suitible for pumping into the
      % FPGA is to build a matrix and then reshape.
      for k=1:self.N_paths
        xt = zeros(self.samps_per_bounce,1);
        yt = xt;
        met_idx = xt*0+k;

        % the setpoint has a meta-idx=0;
        vec_k = [xt(1);
                 yt(1);
                 0];
        
        met_idx(end) = -1;
        vec_k = [vec_k, [xt(:)'; yt(:)'; met_idx(:)']];  %#ok<AGROW>
        vec = [vec; reshape(vec_k, [], 1)];              %#ok<AGROW>
      end
      
    end
  end
end