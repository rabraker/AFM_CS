dat = csvread('C:\Users\arnold\Documents\labview\afm_imaging\data\force-map-04.csv');

%%
volt2nm = (7/20)*1000;
err = dat(:,1);
uz = dat(:,2);

po = err(1)*0;

p = err - po;

figure(1); hold on
plot(p);
ylabel('p')

figure(2);
plot(uz);
ylabel('uz')

figure(3)
plot(uz, p)

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





