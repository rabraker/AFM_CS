Ts = 0.01;

f = 1;
w1 = 2*pi*f

t = (0:1:200)'*Ts;

y = sin(w1*t);
x = y + sin(w1*2*t);

f1 = figure(1); clf
plot(t, x)

save_fig(f1, 'funky_sin')