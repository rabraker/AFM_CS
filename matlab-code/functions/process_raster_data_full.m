function [ output_args ] = process_raster_data_full(dat_path, parent_path, Ts, npix, width )

volts2micron = 50/10;
micron2pix = npix/width;
volts2pix = volts2micron * micron2pix;



datmat = csvread(dat_path);
parent_dat = csvread(parent_path);

xyref = reshape(parent_dat', 2, [])';
xref = xyref(:,1);


samps_per_period = size(parent_dat,1)/2; % twice as many in here for x & y.
samps_per_line = samps_per_period/2;

% datamat data spec:
% [x_s, y_s, err_s, u_s]
datmat = datmat([1:nperiods*samps_per_period], :);

pixmat2 = bin_raster_slow(datmat(:,[1,2,4]), pix, samps_per_period, volts2pix);



end

