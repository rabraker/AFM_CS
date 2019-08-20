classdef RasterTraj
  
  properties
    width_mic;
    raster_freq;
    points_per_line;
    n_lines;    
    x_rate_volts_per_tick;
    x_rate_mic_per_sec;
    x_rate_pix_per_sec;
    
    x_offset_mic;
    y_offset_mic;
    
    y_traj_volts;
    x_traj_volts;
    x_traj_volts_fourier;
    
  end
  
  methods
    function self = RasterTraj(width_mic, raster_freq, n_lines, varargin)
      opts = inputParser();
      opts.addParameter('x_offset_mic', 0)
      opts.addParameter('y_offset_mic', 0)
      opts.addParameter('y_traj_gen', @self.y_traj_gen)
      opts.parse(varargin{:});
      
      y_TG = opts.Results.y_traj_gen;
      self.x_offset_mic = opts.Results.x_offset_mic;
      self.y_offset_mic = opts.Results.y_offset_mic;
      
      self.width_mic = width_mic;
      self.n_lines = n_lines;
      raster_period = 1/raster_freq;
      raster_amplitude = width_mic/2;
      
      square_num_lines = n_lines;
      y_height = (1/square_num_lines)*width_mic;
      
      
      [x_rasterdata, true_freq, points_per_line] = raster(raster_freq, AFM.Ts, 1, 'coerce', 1,...
          'shift', -raster_period/4);
      
      self.raster_freq = true_freq;
      self.points_per_line = points_per_line;

      % we want to start at 0, not -1.
      x_rasterdata.Data = (x_rasterdata.Data+1)*raster_amplitude;
      
      x_rasterdata.Data = x_rasterdata.Data + self.x_offset_mic;
      
      %y_rasterdata = timeseries(linspace(0, y_height, length(x_rasterdata.Time))',...
          %x_rasterdata.Time);
      y_rasterdata = y_TG(y_height, self.points_per_line*2);
      y_rasterdata.Data = y_rasterdata.Data + self.y_offset_mic;
      
      self.x_traj_volts = x_rasterdata;
      self.y_traj_volts = y_rasterdata;
      
      self.y_traj_volts.Data = self.y_traj_volts.Data * AFM.mic2volt_xy;
      self.x_traj_volts.Data = self.x_traj_volts.Data * AFM.mic2volt_xy;
    end
    function y_rasterdata = y_traj_gen(~, y_height, N)
        t = (0:N-1)'*AFM.Ts;
        y = linspace(0, y_height, N)';
        y_rasterdata = timeseries(y, t);
    end
    function [u_ff, u_des] = fourier_truncate(self, n_harmonics, Hyr, w_max)
        
        uu = self.x_traj_volts;
        N = length(uu);
        t = double(0:1:N-1)'*(AFM.Ts);
        
        u_des = 0;
        u_ff = 0;
        for k=1:2:2*n_harmonics
            wk = k*self.raster_freq*2*pi;
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
   
  end
end