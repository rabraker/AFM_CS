
% The problem of connecting the cantilevar end to the sample surface is the
% same is finding where a circle and line intersect. Show that we did this
% correctly.


clc
ell = 3;
yb = 1.5;
b = -0.1275;

surface_width  = 5;
surface_height = 1;

cant_len = 3;
cant_base_ht = 1.5;
tip_ht = 0.4;
tip_wd = 0.3;

m = surface_height/surface_width; % slope

x = -ell:0.001:ell;
% taking this negative part means the cantilever hangs down, not up.
yl = @(x) (2*yb - sqrt( (2*yb)^2 - 4*( x.^2 + yb^2 - ell^2) ) )/2;

figure(2); clf
plot(x, yl(x))
hold on
plot(x, yl(x))
plot(x, m*x +b)
grid 


% 0 = xl^2 [1- m^2] -2*m*[yb - b] * xl + yb^2 -2yb*b + b^2 - l^2
% 0 = xl^2 - 2*m*[yb - b]           yb^2 -2yb*b + b^2 - l^2
% =0        ------------ * xl    + ------------------------
%              [1- m^2]                [1- m^2]

%   xl = -B +- sqrt(B^2 - 4ac)/(2a) . Dont overload b.
yb = cant_base_ht;
ell = cant_len;
a = (1+ m^2); % ok
B = 2*m*( b - yb );
c = yb^2 -2*yb*b + b^2 - ell^2;

xl_1 = ( -B + sqrt( B^2 - 4*a*c) ) / (2 * a);
% this should be negative. Where the line intersects the circle for xl
% negative.
xl_2 = ( -B - sqrt( B^2 - 4*a*c) ) / (2 * a) ;

yy = yl(xl_1);

plot(xl_1, yy, 'x')





F = figure(1); clf
set(F, 'Position', [-1286 389 600 200], 'Color', 'w');
ax = axes('Position', [0.0610 0.1100 0.8913 0.8535]);
ax.NextPlot = 'replaceChildren';

ax2 = axes('Position', [0.5864 0.1400 0.3300 0.4745], 'XTick', [], 'YTick', []);
ax2.NextPlot = 'replaceChildren';


leg = legend(ax2, []);
set(leg, 'NumColumns', 2, 'Position', [0.5630 0.3480 0.1717 0.0370])


figure(F);
F.CurrentAxes = ax;
hold on
xlim([-5, 15])
ylim([-.1, 5])

Ts = 0.01;
t = (0:200)'*Ts;
w = 2*pi;

x = t*0;
z = x;

Frames(length(t)) = struct('cdata',[],'colormap',[]);

for k=1:length(t)
  xo = sin(w*Ts*k);
    
  set(ax2, 'XTick', [], 'YTick', [])
  x_tip = 2.5;
%   y_tip = (x_tip + xo)* tan(theta);
  
%   ytip_base = y_tip + tip_ht;
  
  % The intercept moves as our triangle moves.
  b = -m * xo; % intercept.
  
  % solve quadratic. x = [-b +-sqrt(b^2 - 4*c) ]/2/a. Dont overload b.
  %   xl = -B +- sqrt(B^2 - 4ac)/(2a) . Dont overload b.
  yb = cant_base_ht;
  ell = cant_len;

  a = (1+ m^2); % ok
  B = 2*m*( b - yb );
  c = yb^2 -2*yb*b + b^2 - ell^2;
 
  % Would be more robust to figure out how switch on the negative root.
  xl_1 = ( -B + sqrt( B^2 - 4*a*c) ) / (2 * a);
  % This should be negative. Where the line intersects the circle for xl
  % negative.
%   xl_2 = ( -B - sqrt( B^2 - 4*a*c) ) / (2 * a);
  
  hands1 = draw_right_triangle(ax, xo, surface_width, surface_height);
  hands2 = draw_afm_head(ax, yl(xl_1), xl_1, tip_wd, tip_ht);
  h_cant = draw_cantilevar(ax, yl(xl_1)+tip_ht, xl_1, cant_base_ht);
  % return
  
  ylim(ax, [-0.1, 2.5])
  
  x(k) = xo;
  z(k) = yl(xl_1);
  
  h1 = plot(ax2, t(1:k),  x(1:k), 'b');
  hold(ax2, 'on')
  h2 = plot(ax2, t(1:k), z(1:k)-z(1), 'r');
  xlim(ax2, [0, t(end)])
  ylim(ax2, [-1, 1]*1.05);
  set(h1, 'DisplayName', '$x(t)$')
  set(h2, 'DisplayName', '$z(t)$')
  leg = legend(ax2, [h1, h2]);
  set(leg, 'NumColumns', 2, 'Position', [0.5864 0.6098 0.2566 0.1398],...
    'FontSize', 12)


  drawnow();
  Frames(k) = getframe(F);
  

  if k < length(t)
    delete(hands1)
    delete(hands2)
    delete(h_cant)
    delete(h1);
    delete(h2);
  end
end
%%
% file extension is added automatically. Matlab is shitty and wont do MP4 on
% linux. Convert later with 
% >>> ffmpeg -i sample_tilt.avi sample_tilt.mp4
v = VideoWriter('~/matlab/afm-cs/matlab-code/notes/figures/sample_tilt')


open(v)
for k=1:length(Frames)
  v.writeVideo(Frames(k))
end
v.close()

%%
function h = draw_cantilevar(ax, y_tip_base, x_tip, cant_base_ht);
  
  h = plot(ax, [0, x_tip], [cant_base_ht, y_tip_base], 'k',...
    'LineWidth', 2);
  
  
  
end

function hands = draw_afm_head(ax, y_tip, x_tip, wd, ht)
  
  h1 = plot(ax, [x_tip, x_tip + wd/2], [y_tip, y_tip+ht], 'k');
  h2 = plot(ax, [x_tip, x_tip - wd/2], [y_tip, y_tip+ht], 'k');
  h3 = plot(ax, [x_tip - wd/2, x_tip+wd/2], [y_tip+ht, y_tip+ht], 'k');
  
  hands = [h1, h2, h3];
end

function [hands] = draw_right_triangle(ax, xo, width, height)
  
  h = sqrt(width^2 + height^2);
  
  hbase = plot(ax, [xo, xo+width], [0, 0], '-k');
  hside = plot(ax, [xo, xo]+width, [0, height], '-k');
  hhyp = plot(ax, [xo, xo+width], [0, height], '-k');
  
  hands = [hbase, hside, hhyp];
end