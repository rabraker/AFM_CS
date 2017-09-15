
clear; clc

matroot = getMatPath();
notepath = '/home/arnold/gradschool/proposals/nsf_afm_sept2017/latex/figs';
model_path = fullfile(matroot, 'publications', 'afmMPC', 'mimo_modelData.mat');
load(model_path);

PLANT_preT_x.InputDelay = 0;
G_x = PLANT_preT_x;
Xss  = (eye(size(G_x.B, 1)) - G_x.A)\G_x.B;
Ts = 40e-6;
Ki_x = 0.01;
D_x = ss(tf([1]*Ki_x, [1 -1], Ts));
H = feedback(D_x*G_x, 1);
% --------------------------------------------

mpl = load('upathlocations.txt');

x = mpl(:,1) ;
% - mpl(1,1);
y= mpl(:,2);
% - mpl(1,2);;
lngs = mpl(:,3);
dirs = mpl(:,4);

figure(100)
subplot(2,1,1)
ax1 = gca;
plot(x)

subplot(2,1,2)
plot(diff(x))
ax2 = gca;
linkaxes([ax1, ax2], 'x')

figure(200);  hold on;
plot(x,y)
xlim([0, 512])
ylim([0, 512])
xlm = xlim;
ylm = ylim;

xgrid = repmat([0:1:512]-.5, 2, 1);
ygrid = xgrid;
xg_plus = repmat(xlm(2), 1, 513);
xg_min = repmat(xlm(1), 1, 513);

yg_plus = repmat(ylm(2), 1, 513);
yg_min = repmat(ylm(1), 1, 513);

h1 = plot([xg_min; xg_plus], xgrid, 'Color',  [.82 .82 .78]);

plot(ygrid, [yg_min; yg_plus], 'Color',  [.82 .82 .78]);
ax3 = gca;
xlim([0,100])
ylim([0, 100])

%%
xnew = x;
ynew = y;
k = 1;

while 1
    dx = xnew(k+1) - xnew(k);
    if abs(dx) >1
%        keyboard 
        xnew = [xnew(1:k); xnew(k)+[1:1:dx-1]'; xnew(k+1:end)];
%         ynew = [ynew(1:k); linspace(ynew(k), ynew(k+1), dx-1)'; ynew(k+1:end)];
    end
    
    k = k+1;
    if k == length(x)
        break
    end
end

figure(300); clf; 
plot(diff(xnew))
hold on;
plot(diff(x), '--')

[X, f_s] = power_spec(x, 1);
[Xnew, f_s_new] = power_spec(xnew, 1);

figure(10);
semilogx(f_s, 20*log10(X))
hold on
semilogx(f_s_new, 20*log10(Xnew))


%%
T = 1;
mu2volts = 1/5;
pix2mu = 10/512;
secs_per_pix = (T/2)/512;

t_pix = linspace(0, secs_per_pix*length(x), length(x));
t_mu = [0:Ts:t_pix(end)]';

x_mu = interp1(t_pix, x, t_mu)*pix2mu*mu2volts;


y1 = lsim(H, x_mu, t_mu);
% ---------------------
figure;
plot(t_mu, y1);
%%
t_pix = linspace(0, secs_per_pix*length(xnew), length(xnew));
t_mu = [0:Ts:t_pix(end)]';

x_mu = interp1(t_pix, xnew, t_mu)*pix2mu*mu2volts;


y2 = lsim(H, x_mu, t_mu);

hold on
plot(t_mu, y2-x_mu, '--')


