function x_pixelized = pixelize(x, npix)
    
    kk = (floor(length(x)/npix))*npix;
    
    X = reshape(x(1:kk), [], npix);
    x_pixelized = mean(X)';
    if kk < length(x)
        x_pixelized = [x_pixelized; mean(x(kk+1:end))];
    end

end