function self = load_raw_data(self, raster_paths, opts)
% Loads the raw raster data from the paths defined in raster_paths. This is
% a helper function for the main constructor.
% 
%Arguments
%---------
%  raster_paths :  struct of paths
%                  Generate with the function get_raster_paths(). This 
%                  function expects that cs-paths will minimally contain 
% 
%  raster_paths.meta_path : (string) 
%                           Path to the matfile saved by laview and contains meta
%                           data generated in the experiment, like
%                           control parameters and settings.
% 
%  raster_paths.parent_meta_path : (string)
%                                  This is the .mat file containing the parents
%                                  meta data, saved, e.g., by get_raster_paths.m.
%   
%  raster_paths.data_path : This is the csv file of the actual raw data, saved
%                           by labview.
%  opts :  struct
%          an options structure (typically from name-value pairs
%          passed to the constructor). We exect the following
%          properties.
%  opts.channel_map ChannelMap object.
%                   Defines the indeces for x,y,ze,uz data
%                   channells in the csv file.
%  opts.Ts : double
%            If empty, will use the value in AFM.Ts
%  opts.feature_height: double
%                       Expected feature height of the sample. Used
%                       for visualizing when verbose=true in
%                       process_raster_data().
%  opts.gg :  lti model
%             If this is non-empty, the uz control data will be
%             filtered through this transfer function. Useful for
%             removing drift etc.
% 
%Returns
%-------
%  self : RasterExp
%         The current instance with data loaded into the proper
%         properties.    
% 
%See Also
%--------
%  get_raster_paths.m, RasterExp.load_mat()
% 
  
  channel_map = opts.channel_map;
  self.channel_map = opts.channel_map;
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
  self.yref = xyref(:,2); 
  
  npix = meta_data.Cluster.raster_meta_in.npix;
  width = meta_data.Cluster.raster_meta_in.width;
  self.npix = npix;
  self.width = width;


  % Do a sanity check on width:
  xmax = max(datmat(:, channel_map.x));
  xmin = min(datmat(:, channel_map.x));
  % xref is volts. Convert to microns.
  width_exp = round( (xmax - xmin) * AFM.volts2mic_xy, 0);
  
  if abs(width_exp - width) > .2
    warning(['Experiment width is %f, but claimed width from meta data is %f'....
      'Are you loading the right files?'], width_exp, width);
  end
  
  % Compute unit conversions
  micron2pix = npix/width;
  volts2pix = AFM.volts2mic_xy  * micron2pix;
  self.volts2pix = volts2pix;
  self.micron2pix = micron2pix;
  
  self.samps_per_period = size(parent_dat,1)/2; % twice as many in here for x & y.
  self.samps_per_line = self.samps_per_period/2;
  if floor(self.samps_per_line)~= self.samps_per_line
    warning('Non integer number of samples per line = %f.', ...
            self.samps_per_line)
  end
  
  % Pull out data. First, drop anything extra that got collected.
%   datmat = datmat([1:npix*self.samps_per_period], :);

self.x = datmat([1:npix*self.samps_per_period], channel_map.x);
self.y = datmat([1:npix*self.samps_per_period], channel_map.y);
self.ze = datmat([1:npix*self.samps_per_period], channel_map.ze);
self.uz = datmat([1:npix*self.samps_per_period], channel_map.uz);


  
end
