function rdat = raster_full_to_img(dat_path, parent_path, rdat)

    rdat.datmat = csvread(dat_path);
    rdat.parent_dat = csvread(parent_path);

    xyref = reshape(rdat.parent_dat, 2, [])';
    rdat.xref = xyref(:,1);
%     rdat.Ts = Ts;
%     rdat.npix = npix;
%     rdat.volts2pix = volts2pix;
%     rdat.width = width;
    
    rdat.samps_per_period = size(rdat.parent_dat,1)/2; % twice as many in here for x & y.
    rdat.samps_per_line = rdat.samps_per_period/2;
    rdat.datmat = rdat.datmat([1:rdat.npix*rdat.samps_per_period], :);

    [pixmat, pixelifsampled] = bin_raster_slow(rdat.datmat(:,[1,2,4]),...
                                rdat.npix, rdat.samps_per_period, rdat.volts2pix);
    I_fit = detrend_sampled_plane(pixmat, pixelifsampled);
    I_fit = (I_fit - min(min(I_fit))).*pixelifsampled;
    rdat.pixmat = pixmat;
    rdat.I_fit = I_fit;
    rdat.pixelifsampled = pixelifsampled;
%     rdat.freq = freqs(i);
%     raster_dat_s = [raster_dat_s, rdat];



end