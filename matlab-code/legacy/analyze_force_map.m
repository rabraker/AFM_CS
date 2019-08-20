clear
clc
% dat2 = csvread('C:\Users\arnold\Documents\labview\afm_imaging\data\force-map-04.csv');
addpath('functions');

% dat_s{1} = csvread('C:\Users\arnold\Documents\labview\afm_imaging\data\forcemap_data-out02.csv')
% dat_s{2} = csvread('C:\Users\arnold\Documents\labview\afm_imaging\data\forcemap_data-out03.csv')
% dat_s{2} = csvread(fullfile(getdataroot, 'force-map', 'force-map-03.csv'));

% dat_s{1} = csvread('/media/labserver/afm-cs/force-curve-10-8-2018-04.csv')
% dat_s{1} = csvread('/media/labserver/afm-cs/force-curve-10-8-2018-05.csv')
dat_s{1} = csvread('/media/labserver/afm-cs/force-curve-10-8-2018-06.csv')

Ts = 40e-6;
clf
volt2nm = (7/20)*1000;

figure(3), clf
hold on
figure(4), clf, 
ax2=subplot(2,1,1); 
hold on
ax3=subplot(2,1,2); 
hold on



for datt = dat_s
  figure(3);
    dat = datt{1};
err = dat(:,1);
uz = dat(:,2);

[~, idx_min] = min(uz);

plot(uz, err)

figure(4)
subplot(2,1,1)
t = (0:length(err)-1)'*Ts;
plot(t, err)
subplot(2,1,2)
plot(t, uz)
end

figure(3)
xlabel('uz [v]')
ylabel('p = err - po [v]')
grid on

figure(4)
subplot(2,1,1)
xlabel('t [s]')
ylabel('z-err')
grid on

subplot(2,1,2)
xlabel('t [s]')
ylabel('u_z')
grid on

linkaxes([ax2, ax3], 'x')
%%

uz_up = uz(idx_min:end);
err_up = err(idx_min:end);
t_up = t(idx_min:end);
figure(5), clf

subplot(2,1,1)
plot(t_up, err_up)
xlabel('t [s]')
ylabel('z-err')

grid on, hold on
subplot(2,1,2)
plot(t_up, uz_up)
grid on, hold on
xlabel('t [s]')
ylabel('u-z')


N = 10;

derr_k = 1;
for k=N+1:length(err_up)
  
 derr_k_min1 = derr_k;
 derr_k = err_up(k) - err_up(k-N); 
 
 subplot(2,1,1)
 if derr_k_min1 <0 && derr_k > 0
   plot([t_up(k-N), t_up(k)], [err_up(k-N), err_up(k)], 'r')
   break
 else
   plot([t_up(k-N), t_up(k)], [err_up(k-N), err_up(k)], 'k')
 end
%  if k>77
%   keyboard
%  end
end













%%
k1 = 11599;
k2 = 18546;

% p = [x 1][m b]'

uz_lin = uz(k1:k2)*volt2nm;
p_lin = p(k1:k2);

mb = [uz_lin, 0*uz_lin+1]\p_lin

figure(3); hold on
plot(uz_lin, mb(1)*uz_lin + mb(2));

dfl_volt2nm = abs(1/(mb(1)));

p_nm = p*dfl_volt2nm;

figure(4);
plot(uz*volt2nm, p_nm)

k = mean([.02, .8])*1e-9;

figure(5)
plot(uz*volt2nm, p_nm*k)





