% Fit a line to a set of data and subtract it off.

function x_flat = detrend(x)

    xs = [0:1:length(x)-1]';
    
    PHI = [xs, 0*xs+1];
    
    % x = PHI*[m; b]
    
    mb = PHI\x;
    
    m = mb(1);
    b = mb(2);
    
    x_flat = x - (m*xs + b);

end