clc
dat = csvread('C:\Users\arnold\Documents\labview\afm_imaging\data\data-out2.csv');
Ts = 40e-6;
% [ x(k), y(k), e_z(k),  z(k), u_z(k), FIFO_INDEX(k),  ]
xdat = dat(:,1);
ydat = dat(:,2);
err_dat = dat(:,3);
uz_dat = dat(:,5);
index_dat = dat(:,6);

k = 201;

err_dat = err_dat(k:end);
uz_dat = uz_dat(k:end);
uz_dat = uz_dat - uz_dat(1);
t = [0:1:length(uz_dat)-1]'*Ts;

iddat = iddata(err_dat, uz_dat, Ts);
G = tfest(iddat, 3, 'Ts', Ts)

y = lsim(G, uz_dat, t);

figure(2); clf
subplot(2,1,1)
plot(t, uz_dat)

subplot(2,1,2); hold on
plot(t, err_dat)
plot(t, y)
