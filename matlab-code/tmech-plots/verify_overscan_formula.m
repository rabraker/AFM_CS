addpath functions/state_space_x/

xdirControl = get_xdir_loop_shaped_control();
%%
L = xdirControl.Loop;

Hyrp = xdirControl.Hy_rprime;
Hyr = xdirControl.Hyr;

figure(10),
bode(Hyr)
%%
derz = zpk(1, [], 1, AFM.Ts)

kv = dcgain(minreal((Hyr-1)/derz))

N = floor(1/kv)

v = 4;
t = (0:400)'*AFM.Ts;
u = v*t;
y = lsim(Hyr, u, t);

figure(11); clf
plot(t, y)
hold on
k1 = 200;
plot(t(k1), y(k1), 'x')
plot(t(k1), u(k1), 'x')

u(k1) - y(k1)
ess = v/(kv/AFM.Ts)

idx = find(y > u(k1), 1, 'first') - k1