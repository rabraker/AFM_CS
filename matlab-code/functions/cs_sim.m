

function [bp_im, pixelifsampled] = cs_sim(rdat, csin_data_path)


% Backout the pixelifsampled matrix from the raw csv trajectory data. 
pixelifsampled = cs2pixelifsampled(csin_data_path, rdat.npix, rdat.pix_per_volt);
% Now, restrict only only the data we have (due to imperfect raster
% tracking).
[kmin, kmax] = find_raster_extents(rdat);
pixelifsampled = pixelifsampled(:, kmin:kmax);


I_sampled = rdat.I_fit.*pixelifsampled;


root = '/home/arnold/gradschool/afm-cs';
reconstruct_root = fullfile(root, 'reconstruction/BP');
reconstruct_path = genpath(reconstruct_root);
addpath(reconstruct_path);

% ********* BP *************
bp_im = bp_cs_recon(I_sampled, pixelifsampled);

end