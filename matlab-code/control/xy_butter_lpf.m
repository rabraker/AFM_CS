
fs = 25e3;
wn = 150/(0.5*fs);

[z, p, k] = butter(2, wn, 'low');

gf = zpk(z, p, k, 1/fs);

% gf = gf/dcgain(gf);

figure(1); clf
opts = bodeoptions;
opts.FreqUnits = 'Hz';
bodeplot(gf, opts)
grid
% subplot(2,1,1)
% ylim([-400, 5])

figure(2), pzplot(gf)

figure(3)
step(gf)
%%
s = 200 * 2 * pi;

z = exp(-s * AFM.Ts);

g1 = zpk([0], [z], 1, AFM.Ts)

g = g1*g1*g1;
g = g/dcgain(g);
format shortg
bandwidth(g)/2/pi
% bode(g)
format longg
[num, den] = tfdata(g, 'v')

%%
wz1 = 214 * 2 * pi;
wz2 = 505 * 2 * pi;
zz = 0.01;


g = tf([1, 2*zz*wz1, wz1^2], conv([1, 300*2*pi],  [1, 300*2*pi]));
g2 = tf([1, 2*zz*wz2, wz2^2], conv([1, 450*2*pi],  [1, 450*2*pi]));

g = c2d(g, AFM.Ts, 'matched') * c2d(g2, AFM.Ts, 'matched')
g = tf(g) * (1/dcgain(g))
%%
format shortg
bandwidth(g)/2/pi
% bode(g)
format longg
[num, den] = tfdata(g, 'v')

figure
bode(g)