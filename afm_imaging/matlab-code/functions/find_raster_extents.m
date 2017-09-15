% [kmin, kmax] = find_raster_extents(rdat)
%
% Find the x-axis extent of the raster data. kmin and kmax are the indeces
% where the raster data actually exists. 

function [kmin, kmax] = find_raster_extents(rdat)


    kmin = 1;
    kmax = size(rdat.pixelifsampled, 2);


    for i=1:rdat.npix
       if mean(rdat.pixelifsampled(:,i)) < 0.5
           kmin = i+1;
       else
           break
       end
    end

    for i = kmin+1:rdat.npix
       if mean(rdat.pixelifsampled(:,i)) < 0.5
           kmax = i-1;
           break
       end
    end




end

