

function pixelifsampled = cs2pixelifsampled(cs_data_in_path, npix, pix_per_volt)
cs_data_in = csvread(cs_data_in_path);

cs_dat = reshape(cs_data_in', 3, [])';

pixelifsampled = zeros(npix, npix);
for k = 1:max(cs_dat(:,end))
    % Get indexes corresponding to measurement entity k
   inds = find(cs_dat(:,end) == k); 
   
   % Slice out the data for the current mu-path.
   U_ks = ones(length(inds), 1);
   Y_ks = cs_dat(inds, 2)*pix_per_volt;
   X_ks = cs_dat(inds, 1)*pix_per_volt;
    
   
   y_pix = floor(mean(Y_ks));
   
   x_spread = max(X_ks) - min(X_ks);
   xpix_start = floor(min(X_ks));
   npix_path_k = ceil(x_spread); % number of bins for this path.
   
   % Now, define a set of bins for the x-direction. Each bin will have a
   % different number of data points. 
   % If we include 0, then there are three points for two pixel bins.
   xbins = linspace(min(X_ks), max(X_ks), npix_path_k+1); 
   for jj = 1:npix_path_k
       % Get the indeces correspondiing to the current x-data bin.
       ind_x = find(X_ks >= xbins(jj) & X_ks < xbins(jj+1));
       if ~isempty(ind_x) % Avoid errors if it is empty.
           % Slice out the corresponding height data, and average it
           % together since this is a single pixel.
           u_pix_jj = mean(U_ks(ind_x)); 
           
           % The image index for the y-direction (rows) is y_pix.
           % For the x-direction, the start of the path is at xpix_start.
           % We are inching along the mu-path, so at each iteration of this
           % loop, increment the pixel index by 1, so just use loop iterator,
           % jj.
           I(y_pix+1, xpix_start+jj) = u_pix_jj;
           pixelifsampled(y_pix+1, xpix_start+jj) = 1;
           % Collect the uz data just for plotting/visualization
       end
   end
end

end