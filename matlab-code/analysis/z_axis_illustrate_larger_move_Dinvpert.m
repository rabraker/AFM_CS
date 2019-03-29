% Illustrate that the larger oscillations seen after tip-descent may well just
% be from the larger move, than the 20nm move when scanning over a hole.

num = [1.0  -1.99625,   0.999089];
den = [1.0  -1.99624,   0.999139];
%%
g = tf(den(1:3), num(1:3), AFM.Ts)
damp(g)
figure(300); clf
bode(g)

p = pole(g);
z = zero(g);


w_p = abs(log(p(1)))/AFM.Ts
z_p = -cos(angle(log(p(1))))

w_z = abs(log(z(1)))/AFM.Ts
z_z = -cos(angle(log(z(1))))

zpert = 1 - 0.005;
ppert = 1 + 0.01*0;
g_sp = tf([1, 2*w_z*z_z*zpert*0.95, (w_z*zpert)^2], [1, 2*w_p*z_p*ppert*0.95, (w_p*ppert)^2])

g_zp = c2d(g_sp, AFM.Ts, 'matched')

hold on

bode(g_zp)

%
Ki = 0.1;
pp = 1000*2*pi;
G = tf(pp, [1, pp])
Gz = c2d(G, AFM.Ts, 'matched') * g_zp;
Di = zpk([0], [1], Ki, AFM.Ts);
Dz = Di * (1/g);

L = Dz * Gz;


H = feedback(Di, (1/g)*Gz);
bode(H)

t = [0:1500]'*AFM.Ts;
u = 1 + 0*t;

%%
figure(301); clf

[y, t, x] = lsim(ss(H), u, t);
x0 = (-5)*x(end, :)';
plot(t, y)
hold on

[y, t, x] = lsim(ss(H), u, t, x0);

plot(t, y)
ylim([0, 1.1])

