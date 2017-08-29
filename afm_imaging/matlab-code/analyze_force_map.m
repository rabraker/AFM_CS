
% dat2 = csvread('C:\Users\arnold\Documents\labview\afm_imaging\data\force-map-04.csv');

dat_s{1} = csvread(fullfile(getdataroot, 'data', 'force-map', 'force-map-short_cant-01.csv'));
dat_s{2} = csvread(fullfile(getdataroot, 'data', 'force-map', 'force-map-short_cant-02.csv'));
dat_s{3} = csvread(fullfile(getdataroot, 'data', 'force-map', 'force-map-short_cant-03.csv'));
%%
volt2nm = (7/20)*1000;
figure(3)
hold on
for datt = dat_s
    dat = datt{1};
err = dat(:,1);
uz = dat(:,2);
uz = uz - uz(1);

po = err(1);
p = err - po;

plot(uz, p)


end
% err2 = dat2(:,1);
% uz2 = dat2(:,2);
% uz2 = uz2 - uz2(1);


% po2 = err2(1);

% p2 = err2 - po2;
% figure(1); hold on
% plot(p);
% ylabel('p')
% 
% figure(2);
% plot(uz);
% ylabel('uz')


xlabel('uz [\approx nm]')
ylabel('p = err - po [v]')

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





