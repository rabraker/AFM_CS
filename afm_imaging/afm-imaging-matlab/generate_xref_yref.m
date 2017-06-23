% The purpose of this script is to work out how to translate raster scan
% parameters, e.g. lines/sec, pixel density and scan size into
% ramp rates and triangle wave frequency/amplitude.


clc, clear
% Define parameters:

fs = 25e3;  % sample freq
Ts = 1/fs;


pixel_density = 24;  % always assume square images
% define only a corner, and delta.
x_min = -5;
y_min = -.5;
x_range = 10;

scan_rate = 0.5;  % lines/sec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y_range = x_range; % enforce squareness.

% We only use the trace or re-trace, but not both. 
total_lines = pixel_density*2;


x_lines = pixel_density;

tx = 1/scan_rate;  % sec/line
ty = tx * total_lines;

dxdt = x_range/tx;  % magnitude of rate of change.
dxdk_mag = dxdt*Ts; % magnitude of rate of change per time-step.

dydt = y_range/ty;  % rate of change.
dydk = dydt*Ts      % rate of change per time-step.

t_vec = [0:1:floor(ty/Ts)]*Ts;  % only for plotting
y_vec = zeros(1, floor(ty/Ts)+1);
x_vec = zeros(1, floor(ty/Ts)+1);

y_vec(1) = y_min;
x_vec(1) = x_min;

dxdk = dxdk_mag;

for k = 1:length(y_vec)-1
    % dxkd changes sign at vertices. Detect when x(k+1) leaves the "box"
    % and then change sign of dxdk.
    xk1 = x_vec(k) + dxdk;
    if xk1 > x_min + x_range || xk1 < x_min
        dxdk = -dxdk;
        xk1 = x_vec(k) + dxdk;
    end
    x_vec(k+1) = xk1;
    
   % dydk never changes. 
   y_vec(k+1) = y_vec(k) + dydk;
end

figure(1)
plot(t_vec, y_vec, '--')

figure(2)
plot(t_vec, x_vec)

figure(3)
plot(x_vec, y_vec)

