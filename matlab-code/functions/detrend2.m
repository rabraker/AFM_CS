% Fit a line to a set of data and subtract it off.

function x_flat = detrend2(x)

    xs = [0:1:length(x)-1]';
    
    PHI = [xs.^2, xs, 0*xs+1];
    
    % x = PHI*[m; b]
    
    mb = PHI\x;
    

    c = mb(1);
    m = mb(2);
    b = mb(3);
%     keyboard
    x_flat = x - (c*xs.^2 + m*xs + b);

end